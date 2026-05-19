#ifndef PRECICE_NO_MPI

#include <Eigen/Core>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "com/SharedPointer.hpp"
#include "com/SocketCommunicationFactory.hpp"
#include "m2n/CoRTCommunication.hpp"
#include "m2n/M2N.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/IntraComm.hpp"
#include "xml/ConfigParser.hpp"
#include "xml/XMLTag.hpp"

using precice::testing::TestContext;

using namespace precice;
using namespace m2n;

BOOST_AUTO_TEST_SUITE(M2NTests)
BOOST_AUTO_TEST_SUITE(CoRT)

namespace {

std::vector<double> expandComponents(const std::vector<double> &values, int valueDimension)
{
  std::vector<double> expanded;
  expanded.reserve(values.size() * valueDimension);
  for (double value : values) {
    for (int dim = 0; dim < valueDimension; ++dim) {
      expanded.push_back(value);
    }
  }
  return expanded;
}

void process(std::vector<double> &data)
{
  for (double &value : data) {
    value += utils::IntraComm::getRank() + 1;
  }
}

void runRoundtrip(const TestContext  &context,
                  CoRTCommunication  &communication,
                  std::vector<double> data,
                  std::vector<double> expectedData,
                  int                 valueDimension)
{
  if (context.isNamed("A")) {
    communication.send(data, valueDimension);
    communication.receive(data, valueDimension);
    BOOST_TEST(testing::equals(data, expectedData));
  } else {
    BOOST_TEST(context.isNamed("B"));
    communication.receive(data, valueDimension);
    BOOST_TEST(testing::equals(data, expectedData));
    process(data);
    communication.send(data, valueDimension);
  }
}

void runCoRTComTest1(const TestContext &context)
{
  BOOST_TEST(context.hasSize(2));

  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 2, testing::nextMeshID()));

  CoRTCommunication communication(std::make_shared<com::SocketCommunicationFactory>(), mesh);

  std::vector<double> scalarData;
  std::vector<double> scalarExpectedData;

  if (context.isNamed("A")) {
    if (context.isPrimary()) {
      mesh->setGlobalNumberOfVertices(10);
      mesh->setVertexDistribution({{0, {0, 1, 3, 5, 7}}, {1, {1, 2, 4, 5, 6}}});

      scalarData         = {10, 20, 40, 60, 80};
      scalarExpectedData = {10 + 2, 4 * 20 + 3, 40 + 2, 4 * 60 + 3, 80 + 2};
    } else {
      scalarData         = {20, 30, 50, 60, 70};
      scalarExpectedData = {4 * 20 + 3, 30 + 1, 50 + 2, 4 * 60 + 3, 70 + 1};
    }
  } else {
    BOOST_TEST(context.isNamed("B"));
    if (context.isPrimary()) {
      mesh->setGlobalNumberOfVertices(10);
      mesh->setVertexDistribution({{0, {1, 2, 5, 6}}, {1, {0, 1, 3, 4, 5, 7}}});

      scalarData.assign(4, -1);
      scalarExpectedData = {2 * 20, 30, 2 * 60, 70};
    } else {
      scalarData.assign(6, -1);
      scalarExpectedData = {10, 2 * 20, 40, 50, 2 * 60, 80};
    }
  }

  if (context.isNamed("A")) {
    communication.requestConnection("B", "A");
  } else {
    communication.acceptConnection("B", "A");
  }

  runRoundtrip(context, communication, scalarData, scalarExpectedData, 1);
  runRoundtrip(context, communication, expandComponents(scalarData, 2), expandComponents(scalarExpectedData, 2), 2);
}

void runNoOverlapTest(const TestContext &context)
{
  BOOST_TEST(context.hasSize(2));

  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 2, testing::nextMeshID()));

  CoRTCommunication communication(std::make_shared<com::SocketCommunicationFactory>(), mesh);

  std::vector<double> data = {static_cast<double>(context.rank + 1)};

  if (context.isNamed("A")) {
    if (context.isPrimary()) {
      mesh->setGlobalNumberOfVertices(4);
      mesh->setVertexDistribution({{0, {0}}, {1, {2}}});
    }
    communication.requestConnection("B", "A");
    communication.send(data, 1);
    communication.receive(data, 1);
  } else {
    BOOST_TEST(context.isNamed("B"));
    if (context.isPrimary()) {
      mesh->setGlobalNumberOfVertices(4);
      mesh->setVertexDistribution({{0, {1}}, {1, {3}}});
    }
    communication.acceptConnection("B", "A");
    communication.receive(data, 1);
    communication.send(data, 1);
  }

  BOOST_TEST(testing::equals(data, std::vector<double>{0.0}));
}

void runSameConnectionTest(const TestContext &context)
{
  BOOST_TEST(context.hasSize(2));

  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 2, testing::nextMeshID()));

  if (context.isPrimary()) {
    mesh->setConnectedRanks({0});
  } else {
    mesh->setConnectedRanks({1});
  }

  CoRTCommunication communication(std::make_shared<com::SocketCommunicationFactory>(), mesh);
  std::vector<int>  receiveData;

  if (context.isNamed("A")) {
    communication.requestPreConnection("Solid", "Fluid");
    communication.broadcastSend(context.isPrimary() ? 5 : 10);
  } else {
    BOOST_TEST(context.isNamed("B"));
    communication.acceptPreConnection("Solid", "Fluid");
    communication.broadcastReceiveAll(receiveData);
    BOOST_TEST_REQUIRE(receiveData.size() == 1);
    BOOST_TEST(receiveData.at(0) == (context.isPrimary() ? 5 : 10));
  }
}

void runCrossConnectionTest(const TestContext &context)
{
  BOOST_TEST(context.hasSize(2));

  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 2, testing::nextMeshID()));

  if (context.isPrimary()) {
    mesh->setConnectedRanks({1});
  } else {
    mesh->setConnectedRanks({0});
  }

  CoRTCommunication communication(std::make_shared<com::SocketCommunicationFactory>(), mesh);
  std::vector<int>  receiveData;

  if (context.isNamed("A")) {
    communication.requestPreConnection("Solid", "Fluid");
    communication.broadcastSend(context.isPrimary() ? 5 : 10);
  } else {
    BOOST_TEST(context.isNamed("B"));
    communication.acceptPreConnection("Solid", "Fluid");
    communication.broadcastReceiveAll(receiveData);
    BOOST_TEST_REQUIRE(receiveData.size() == 1);
    BOOST_TEST(receiveData.at(0) == (context.isPrimary() ? 10 : 5));
  }
}

void runEmptyConnectionTest(const TestContext &context)
{
  BOOST_TEST(context.hasSize(2));

  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 2, testing::nextMeshID()));

  if (context.isPrimary()) {
    mesh->setConnectedRanks({0});
  }

  CoRTCommunication communication(std::make_shared<com::SocketCommunicationFactory>(), mesh);
  std::vector<int>  receiveData;

  if (context.isNamed("A")) {
    communication.requestPreConnection("Solid", "Fluid");
    communication.broadcastSend(5);
  } else {
    BOOST_TEST(context.isNamed("B"));
    communication.acceptPreConnection("Solid", "Fluid");
    communication.broadcastReceiveAll(receiveData);

    if (context.isPrimary()) {
      BOOST_TEST_REQUIRE(receiveData.size() == 1);
      BOOST_TEST(receiveData.at(0) == 5);
    } else {
      BOOST_TEST(receiveData.empty());
    }
  }
}

void runMeshBroadcastTest(const TestContext &context)
{
  BOOST_TEST(context.hasSize(2));

  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 2, testing::nextMeshID()));

  if (context.isNamed("A")) {
    Eigen::VectorXd position(2);
    if (context.isPrimary()) {
      position << 5.5, 0.0;
      mesh::Vertex &v1 = mesh->createVertex(position);
      position << 1.0, 2.0;
      mesh::Vertex &v2 = mesh->createVertex(position);
      mesh->createEdge(v1, v2);
      mesh->setConnectedRanks({0});
    } else {
      position << 1.5, 0.0;
      mesh::Vertex &v1 = mesh->createVertex(position);
      position << 1.5, 2.0;
      mesh::Vertex &v2 = mesh->createVertex(position);
      mesh->createEdge(v1, v2);
      mesh->setConnectedRanks({1});
    }
  } else {
    BOOST_TEST(context.isNamed("B"));
    mesh->setConnectedRanks({context.isPrimary() ? 0 : 1});
  }

  CoRTCommunication communication(std::make_shared<com::SocketCommunicationFactory>(), mesh);

  if (context.isNamed("A")) {
    communication.requestPreConnection("Solid", "Fluid");
    communication.broadcastSendMesh();
  } else {
    communication.acceptPreConnection("Solid", "Fluid");
    communication.broadcastReceiveAllMesh();

    BOOST_TEST(mesh->nVertices() == 2);
    if (context.isPrimary()) {
      BOOST_TEST(mesh->vertex(0).coord(0) == 5.50);
      BOOST_TEST(mesh->vertex(0).coord(1) == 0.0);
      BOOST_TEST(mesh->vertex(1).coord(0) == 1.0);
      BOOST_TEST(mesh->vertex(1).coord(1) == 2.0);
    } else {
      BOOST_TEST(mesh->vertex(0).coord(0) == 1.50);
      BOOST_TEST(mesh->vertex(0).coord(1) == 0.0);
      BOOST_TEST(mesh->vertex(1).coord(0) == 1.50);
      BOOST_TEST(mesh->vertex(1).coord(1) == 2.0);
    }
  }
}

void runCommunicationMapTest(const TestContext &context)
{
  BOOST_TEST(context.hasSize(2));

  mesh::PtrMesh                   mesh(new mesh::Mesh("Mesh", 2, testing::nextMeshID()));
  std::map<int, std::vector<int>> localCommunicationMap;

  if (context.isNamed("A")) {
    if (context.isPrimary()) {
      mesh->setConnectedRanks({0});
      localCommunicationMap[0] = {102, 1022, 10222};
      localCommunicationMap[1] = {103, 1033, 10333};
    } else {
      mesh->setConnectedRanks({1});
      localCommunicationMap[0] = {112, 1122, 11222};
      localCommunicationMap[1] = {113, 1133, 11333};
    }
  } else {
    BOOST_TEST(context.isNamed("B"));
    mesh->setConnectedRanks({context.isPrimary() ? 0 : 1});
  }

  CoRTCommunication communication(std::make_shared<com::SocketCommunicationFactory>(), mesh);

  if (context.isNamed("A")) {
    communication.requestPreConnection("Solid", "Fluid");
    communication.scatterAllCommunicationMap(localCommunicationMap);
  } else {
    communication.acceptPreConnection("Solid", "Fluid");
    communication.gatherAllCommunicationMap(localCommunicationMap);

    BOOST_TEST(localCommunicationMap.size() == 1);
    if (context.isPrimary()) {
      BOOST_TEST_REQUIRE(localCommunicationMap.at(0).size() == 3);
      BOOST_TEST(localCommunicationMap.at(0).at(0) == 102);
      BOOST_TEST(localCommunicationMap.at(0).at(1) == 1022);
      BOOST_TEST(localCommunicationMap.at(0).at(2) == 10222);
    } else {
      BOOST_TEST_REQUIRE(localCommunicationMap.at(1).size() == 3);
      BOOST_TEST(localCommunicationMap.at(1).at(0) == 113);
      BOOST_TEST(localCommunicationMap.at(1).at(1) == 1133);
      BOOST_TEST(localCommunicationMap.at(1).at(2) == 11333);
    }
  }
}

std::map<int, std::vector<int>> localMapForTwoLevelData(const TestContext &context)
{
  if (context.isNamed("B")) {
    if (context.isPrimary()) {
      return {{0, {0, 2}}, {1, {0, 1, 2, 3}}};
    }
    return {{0, {0, 1, 2, 4, 5}}, {1, {1, 3, 4}}};
  }
  return {};
}

std::map<int, std::vector<int>> remoteMapForTwoLevelData(const TestContext &context)
{
  if (context.isPrimary()) {
    return {{0, {1, 3}}, {1, {0, 1, 3, 4}}};
  }
  return {{0, {0, 1, 2, 3, 4}}, {1, {0, 2, 3}}};
}

void runTwoLevelDataExchangeTest(const TestContext &context)
{
  BOOST_TEST(context.hasSize(2));

  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 2, testing::nextMeshID()));
  mesh->setConnectedRanks({0, 1});

  CoRTCommunication communication(std::make_shared<com::SocketCommunicationFactory>(), mesh);

  std::vector<double> scalarData;
  std::vector<double> scalarExpectedData;

  if (context.isNamed("A")) {
    if (context.isPrimary()) {
      scalarData         = {10, 20, 40, 60, 80};
      scalarExpectedData = {10 + 2, 4 * 20 + 3, 40 + 2, 4 * 60 + 3, 80 + 2};
    } else {
      scalarData         = {20, 30, 50, 60, 70};
      scalarExpectedData = {4 * 20 + 3, 30 + 1, 50 + 2, 4 * 60 + 3, 70 + 1};
    }

    communication.requestPreConnection("B", "A");
    mesh::Mesh::CommunicationMap gatheredCommunicationMap;
    communication.gatherAllCommunicationMap(gatheredCommunicationMap);
    mesh->getCommunicationMap() = std::move(gatheredCommunicationMap);
  } else {
    BOOST_TEST(context.isNamed("B"));
    if (context.isPrimary()) {
      scalarData.assign(4, -1);
      scalarExpectedData = {2 * 20, 30, 2 * 60, 70};
    } else {
      scalarData.assign(6, -1);
      scalarExpectedData = {10, 2 * 20, 40, 50, 2 * 60, 80};
    }

    communication.acceptPreConnection("B", "A");
    mesh::Mesh::CommunicationMap remoteCommunicationMap = remoteMapForTwoLevelData(context);
    communication.scatterAllCommunicationMap(remoteCommunicationMap);
    mesh->getCommunicationMap() = localMapForTwoLevelData(context);
  }

  communication.completeSecondaryRanksConnection();

  runRoundtrip(context, communication, scalarData, scalarExpectedData, 1);
  runRoundtrip(context, communication, expandComponents(scalarData, 2), expandComponents(scalarExpectedData, 2), 2);
}

m2n::PtrM2N configureM2N(const std::string &filename)
{
  xml::XMLTag                          root = xml::getRootTag();
  m2n::M2NConfiguration::SharedPointer m2nConfig(new m2n::M2NConfiguration(root));
  const xml::ConfigurationContext      ccontext{"A", 0, 2};
  xml::configure(root, ccontext, testing::getPathToSources() + "/m2n/tests/" + filename);
  return m2nConfig->getM2N("A", "B");
}

} // namespace

PRECICE_TEST_SETUP("A"_on(2_ranks).setupIntraComm(), "B"_on(2_ranks).setupIntraComm(), Require::Events)
BOOST_AUTO_TEST_CASE(SendReceiveRebuildsPatterns)
{
  PRECICE_TEST();
  runCoRTComTest1(context);
}

PRECICE_TEST_SETUP("A"_on(2_ranks).setupIntraComm(), "B"_on(2_ranks).setupIntraComm(), Require::Events)
BOOST_AUTO_TEST_CASE(NoOverlapDoesNotDeadlock)
{
  PRECICE_TEST();
  runNoOverlapTest(context);
}

PRECICE_TEST_SETUP("A"_on(2_ranks).setupIntraComm(), "B"_on(2_ranks).setupIntraComm(), Require::Events)
BOOST_AUTO_TEST_CASE(TwoLevelSameConnection)
{
  PRECICE_TEST();
  runSameConnectionTest(context);
}

PRECICE_TEST_SETUP("A"_on(2_ranks).setupIntraComm(), "B"_on(2_ranks).setupIntraComm(), Require::Events)
BOOST_AUTO_TEST_CASE(TwoLevelCrossConnection)
{
  PRECICE_TEST();
  runCrossConnectionTest(context);
}

PRECICE_TEST_SETUP("A"_on(2_ranks).setupIntraComm(), "B"_on(2_ranks).setupIntraComm(), Require::Events)
BOOST_AUTO_TEST_CASE(TwoLevelEmptyConnection)
{
  PRECICE_TEST();
  runEmptyConnectionTest(context);
}

PRECICE_TEST_SETUP("A"_on(2_ranks).setupIntraComm(), "B"_on(2_ranks).setupIntraComm(), Require::Events)
BOOST_AUTO_TEST_CASE(TwoLevelMeshBroadcast)
{
  PRECICE_TEST();
  runMeshBroadcastTest(context);
}

PRECICE_TEST_SETUP("A"_on(2_ranks).setupIntraComm(), "B"_on(2_ranks).setupIntraComm(), Require::Events)
BOOST_AUTO_TEST_CASE(TwoLevelCommunicationMap)
{
  PRECICE_TEST();
  runCommunicationMapTest(context);
}

PRECICE_TEST_SETUP("A"_on(2_ranks).setupIntraComm(), "B"_on(2_ranks).setupIntraComm(), Require::Events)
BOOST_AUTO_TEST_CASE(TwoLevelDataExchange)
{
  PRECICE_TEST();
  runTwoLevelDataExchangeTest(context);
}

PRECICE_TEST_SETUP("A"_on(1_rank))
BOOST_AUTO_TEST_CASE(ConfigurationAcceptsCoRT)
{
  PRECICE_TEST();
  BOOST_TEST(configureM2N("cort-valid.xml")->usesCoRTCommunication());
}

PRECICE_TEST_SETUP("A"_on(1_rank))
BOOST_AUTO_TEST_CASE(ConfigurationAcceptsTwoLevelInit)
{
  PRECICE_TEST();
  const m2n::PtrM2N configuredM2N = configureM2N("cort-valid-two-level.xml");
  BOOST_TEST(configuredM2N->usesCoRTCommunication());
  BOOST_TEST(configuredM2N->usesTwoLevelInitialization());
}

PRECICE_TEST_SETUP("A"_on(1_rank))
BOOST_AUTO_TEST_CASE(ConfigurationRejectsGatherScatter)
{
  PRECICE_TEST();
  BOOST_CHECK_THROW(configureM2N("cort-invalid-gather-scatter.xml"), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END() // CoRT
BOOST_AUTO_TEST_SUITE_END() // M2NTests

#endif // PRECICE_NO_MPI
