// Some test cases graphs taken from Part 2

#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include "graph.hpp"
#include "my_integer.hpp"

// You can generate some random graphs to help in your testing
// The graph has N vertices and p is the probability there is an
// edge between any two vertices. 
// You can vary seed to get different graphs
Graph<int> randomGraph(int N, unsigned seed, double p = 0.5) {
  std::mt19937 mt {seed};
  // set up random number generator that is 1 with probability p and
  // 0 with probability 1-p
  std::binomial_distribution<int> heads {1, p};
  // Set the minimum and maximum edge weight here
  const int minWeight {1};
  const int maxWeight {10};
  std::uniform_int_distribution<int> weight {minWeight, maxWeight};
  Graph<int> G {N};
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      if (heads(mt)) {
        // add edge between i and j with random weight 
        // between minWeight and maxWeight
        G.addEdge(i, j, weight(mt));
      }
    }
  }
  return G;
}

// You can generate some random graphs to help in your testing
// The graph has N vertices and p is the probability there is an
// edge between any two vertices. 
// You can vary seed to get different graphs
Graph<MyInteger> randomGraphMyInteger(int N, unsigned seed, double p = 0.5) {
  std::mt19937 mt {seed};
  // set up random number generator that is 1 with probability p and
  // 0 with probability 1-p
  std::binomial_distribution<int> heads {1, p};
  // Set the minimum and maximum edge weight here
  const int minWeight {1};
  const int maxWeight {10};
  std::uniform_int_distribution<int> weight {minWeight, maxWeight};
  Graph<MyInteger> G {N};
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      if (heads(mt)) {
        // add edge between i and j with random weight 
        // between minWeight and maxWeight
        G.addEdge(i, j, MyInteger(weight(mt)));
      }
    }
  }
  return G;
}

// singleSourceIndex tests /////////////////////////////////////////////
// Test on int
TEST(DijkstraTest, distanceFrom0IndexTiny) {
  Graph<int> G {"tinyEWD.txt"};
  Graph<int> result = singleSourceIndex<int>(G, 0);
  bool expected = correctShortestPathTree(G, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexMedium) {
  Graph<int> G {"mediumEWD.txt"};
  Graph<int> result = singleSourceIndex<int>(G, 0);
  bool expected = correctShortestPathTree(G, result, 0);
  EXPECT_EQ(expected, true);
}

// Test on double
TEST(DijkstraTest, distanceFrom0IndexTinyDouble) {
  Graph<double> G {"tinyEWD.txt"};
  Graph<double> result = singleSourceIndex<double>(G, 0);
  bool expected = correctShortestPathTree(G, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexMediumDouble) {
  Graph<double> G {"mediumEWD.txt"};
  Graph<double> result = singleSourceIndex<double>(G, 0);
  bool expected = correctShortestPathTree(G, result, 0);
  EXPECT_EQ(expected, true);
}

// Test on MyInteger
TEST(DijkstraTest, distanceFrom0IndexTinyMyInt) {
  Graph<MyInteger> G {"tinyEWD.txt"};
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(G, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(G, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexMediumMyInt) {
  Graph<MyInteger> G {"mediumEWD.txt"};
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(G, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(G, result, 0);
  EXPECT_EQ(expected, true);
}

// Test from different sources
TEST(DijkstraTest, distanceFrom3IndexTinyMyInt) {
  Graph<MyInteger> G {"tinyEWD.txt"};
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(G, 3);
  bool expected = correctShortestPathTree(G, result, 3);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom3IndexMediumMyInt) {
  Graph<MyInteger> G {"mediumEWD.txt"};
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(G, 3);
  bool expected = correctShortestPathTree(G, result, 3);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom6IndexTinyMyInt) {
  Graph<MyInteger> G {"tinyEWD.txt"};
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(G, 6);
  bool expected = correctShortestPathTree(G, result, 6);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom8IndexMediumMyInt) {
  Graph<MyInteger> G {"mediumEWD.txt"};
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(G, 8);
  bool expected = correctShortestPathTree(G, result, 8);
  EXPECT_EQ(expected, true);
}

// Test on USA road
TEST(DijkstraTest, distanceFrom0IndexUSAMyInt) {
  Graph<MyInteger> G {"USA-road-d.NY.gr"};
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(G, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(G, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom10IndexUSAMyInt) {
  Graph<MyInteger> G {"USA-road-d.NY.gr"};
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(G, 10);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(G, result, 10);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom100IndexUSAMyInt) {
  Graph<MyInteger> G {"USA-road-d.NY.gr"};
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(G, 100);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(G, result, 100);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom1000IndexUSAMyInt) {
  Graph<MyInteger> G {"USA-road-d.NY.gr"};
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(G, 1000);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(G, result, 1000);
  EXPECT_EQ(expected, true);
}

// Test on smaller graphs with specific error prone characteristics
TEST(DijkstraTest, distanceFrom0IndexG1) {
  Graph<int> G1 {3};
  G1.addEdge(0, 1, 2);
  G1.addEdge(1, 0, 1);
  G1.addEdge(1, 2, 10);
  G1.addEdge(2, 1, 5);
  Graph<int> result = singleSourceIndex<int>(G1, 0);
  bool expected = correctShortestPathTree(G1, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexG2) {
  Graph<double> G2 {3};
  G2.addEdge(0, 1, 2.1);
  G2.addEdge(1, 0, 1.4);
  G2.addEdge(1, 2, 10.0);
  G2.addEdge(2, 1, 5.89);
  Graph<double> result = singleSourceIndex<double>(G2, 0);
  bool expected = correctShortestPathTree(G2, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexG3) {
  Graph<double> G3 {3};
  G3.addEdge(0, 1, 2.1);
  G3.addEdge(1, 0, 1.4);
  G3.addEdge(1, 2, 10.0);
  G3.addEdge(2, 1, 5.89);
  Graph<double> result = singleSourceIndex<double>(G3, 0);
  bool expected = correctShortestPathTree(G3, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexG4) {
  Graph<double> G4 {2};
  G4.addEdge(0, 1, 2.1);
  G4.addEdge(1, 0, 1.4);
  Graph<double> result = singleSourceIndex<double>(G4, 0);
  bool expected = correctShortestPathTree(G4, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexG5) {
  Graph<double> G5 {3};
  G5.addEdge(0, 1, 2.1);
  G5.addEdge(1, 0, 1.4);
  G5.addEdge(1, 2, 10.0);
  G5.addEdge(2, 1, 5.89);
  Graph<double> result = singleSourceIndex<double>(G5, 0);
  bool expected = correctShortestPathTree(G5, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexG6) {
  Graph<double> G6 {4};
  G6.addEdge(0, 1, 2.1);
  G6.addEdge(1, 0, 1.4);
  G6.addEdge(1, 2, 10.0);
  G6.addEdge(2, 1, 5.89);
  G6.addEdge(3, 1, 1.27);
  Graph<double> result = singleSourceIndex<double>(G6, 0);
  bool expected = correctShortestPathTree(G6, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexG7) {
  Graph<double> G7 {5};
  G7.addEdge(0, 2, 2.1);
  G7.addEdge(0, 4, 1.4);
  G7.addEdge(4, 1, 10.5);
  G7.addEdge(1, 3, 2.7);
  Graph<double> result = singleSourceIndex<double>(G7, 0);
  bool expected = correctShortestPathTree(G7, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexG8) {
  Graph<double> G8 {6};
  G8.addEdge(0, 2, 2.1);
  G8.addEdge(0, 5, 1.4);
  G8.addEdge(5, 1, 10.5);
  G8.addEdge(1, 3, 2.7);
  Graph<double> result = singleSourceIndex<double>(G8, 0);
  bool expected = correctShortestPathTree(G8, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexG9) {
  Graph<double> G9 {10};
  G9.addEdge(0, 2, 2.1);
  G9.addEdge(0, 4, 1.4);
  G9.addEdge(4, 1, 10.5);
  G9.addEdge(1, 3, 2.7);
  G9.addEdge(0, 9, 2.4);
  G9.addEdge(2, 7, 8.8);
  G9.addEdge(5, 6, 9.3);
  Graph<double> result = singleSourceIndex<double>(G9, 0);
  bool expected = correctShortestPathTree(G9, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexG10) {
  Graph<double> G10 {7};
  G10.addEdge(0, 2, 2.1);
  G10.addEdge(0, 4, 1.4);
  G10.addEdge(1, 4, 10.5);
  G10.addEdge(1, 3, 2.7);
  G10.addEdge(0, 5, 2.4);
  G10.addEdge(2, 6, 8.8);
  Graph<double> result = singleSourceIndex<double>(G10, 0);
  bool expected = correctShortestPathTree(G10, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexG11) {
  Graph<double> G11 {7};
  G11.addEdge(0, 2, 2.1);
  G11.addEdge(0, 4, 1.4);
  G11.addEdge(4, 1, 10.5);
  G11.addEdge(1, 3, 2.7);
  G11.addEdge(0, 5, 2.4);
  G11.addEdge(2, 6, 8.8);
  G11.addEdge(6, 5, 6.2);
  Graph<double> result = singleSourceIndex<double>(G11, 0);
  bool expected = correctShortestPathTree(G11, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexG12) {
  Graph<int> G12 {4};
  G12.addEdge(3, 1, 4);
  G12.addEdge(3, 0, 1);
  G12.addEdge(3, 2, 2);
  G12.addEdge(1, 0, 9);
  Graph<int> result = singleSourceIndex<int>(G12, 3);
  bool expected = correctShortestPathTree(G12, result, 3);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexG13) {
  Graph<int> G13 {7};
  G13.addEdge(0, 1, 2);
  G13.addEdge(0, 2, 1);
  G13.addEdge(1, 3, 1);
  G13.addEdge(1, 4, 5);
  G13.addEdge(2, 5, 3);
  G13.addEdge(2, 6, 2);
  Graph<int> result = singleSourceIndex<int>(G13, 0);
  bool expected = correctShortestPathTree(G13, result, 0);
  EXPECT_EQ(expected, true);
}


TEST(DijkstraTest, distanceFrom0IndexG14) {
  Graph<int> G14 {7};
  G14.addEdge(0, 1, 1);
  G14.addEdge(0, 2, 1);
  G14.addEdge(1, 3, 1);
  G14.addEdge(1, 4, 1);
  G14.addEdge(2, 5, 1);
  G14.addEdge(2, 6, 1);
  Graph<int> result = singleSourceIndex<int>(G14, 0);
  bool expected = correctShortestPathTree(G14, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom3IndexG15) {
  Graph<int> G15 {8};
  G15.addEdge(3, 2, 1);
  G15.addEdge(3, 4, 4);
  G15.addEdge(7, 5, 10);
  G15.addEdge(7, 6, 2);
  G15.addEdge(3, 7, 3);
  G15.addEdge(4, 1, 5);
  G15.addEdge(2, 0, 11);
  Graph<int> result = singleSourceIndex<int>(G15, 3);
  bool expected = correctShortestPathTree(G15, result, 3);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexG16) {
  Graph<int> G16 {8};
  G16.addEdge(3, 2, 4);
  G16.addEdge(3, 4, 2);
  G16.addEdge(7, 5, 10);
  G16.addEdge(7, 6, 1);
  G16.addEdge(3, 7, 9);
  G16.addEdge(4, 1, 5);
  G16.addEdge(2, 0, 2);
  Graph<int> result = singleSourceIndex<int>(G16, 0);
  bool expected = correctShortestPathTree(G16, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexG17) {
  Graph<int> G17(1);
  Graph<int> result = singleSourceIndex<int>(G17, 0);
  bool expected = correctShortestPathTree(G17, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexG18) {
  Graph<int> G18(6);
  G18.addEdge(0, 1, 1);
  G18.addEdge(1, 2, 1);
  G18.addEdge(2, 3, 1);
  G18.addEdge(3, 4, 1);
  G18.addEdge(4, 5, 1);
  Graph<int> result = singleSourceIndex<int>(G18, 0);
  bool expected = correctShortestPathTree(G18, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexG19) {
  Graph<int> G19(6);
  G19.addEdge(0, 1, 1);
  G19.addEdge(1, 2, 1);
  G19.addEdge(2, 3, 1);
  G19.addEdge(3, 4, 1);
  G19.addEdge(4, 5, 1);
  G19.addEdge(5, 0, 1);
  G19.addEdge(0, 3, 1);
  Graph<int> result = singleSourceIndex<int>(G19, 0);
  bool expected = correctShortestPathTree(G19, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexG20) {
  Graph<int> G20(4);
  G20.addEdge(0, 1, 4); 
  G20.addEdge(0, 2, 5); 
  G20.addEdge(2, 1, 3); 
  G20.addEdge(1, 3, 1); 
  Graph<int> result = singleSourceIndex<int>(G20, 0);
  bool expected = correctShortestPathTree(G20, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexG21) {
  Graph<int> G21(6);
  G21.addEdge(0, 1, 1); 
  G21.addEdge(0, 2, 1); 
  G21.addEdge(0, 3, 1); 
  G21.addEdge(0, 4, 1); 
  G21.addEdge(0, 5, 1); 
  G21.addEdge(2, 4, 1); 
  G21.addEdge(4, 5, 1); 
  Graph<int> result = singleSourceIndex<int>(G21, 0);
  bool expected = correctShortestPathTree(G21, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexG22) {
  Graph<int> G22(5);
  G22.addEdge(0, 1, 4);
  G22.addEdge(1, 2, 2);
  G22.addEdge(2, 3, 2);
  G22.addEdge(3, 4, 1);
  G22.addEdge(0, 3, 2);
  Graph<int> result = singleSourceIndex<int>(G22, 0);
  bool expected = correctShortestPathTree(G22, result, 0);
  EXPECT_EQ(expected, true);
}

// Test on random graphs of varying sizes
TEST(DijkstraTest, distanceFrom0IndexRandom1) {
  Graph<MyInteger> GR = randomGraphMyInteger(1, 3431);
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexRandom2) {
  Graph<MyInteger> GR = randomGraphMyInteger(2, 3431);
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexRandom4) {
  Graph<MyInteger> GR = randomGraphMyInteger(4, 3431);
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexRandom8) {
  Graph<MyInteger> GR = randomGraphMyInteger(8, 3431);
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexRandom16) {
  Graph<MyInteger> GR = randomGraphMyInteger(16, 3431);
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexRandom32) {
  Graph<MyInteger> GR = randomGraphMyInteger(32, 3431);
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexRandom64) {
  Graph<MyInteger> GR = randomGraphMyInteger(64, 3431);
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexRandom128) {
  Graph<MyInteger> GR = randomGraphMyInteger(128, 3431);
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

// Perform tests on same sized random graphs with less edges
TEST(DijkstraTest, distanceFrom0IndexRandom1Less) {
  Graph<MyInteger> GR = randomGraphMyInteger(1, 3431, 0.2);
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexRandom2Less) {
  Graph<MyInteger> GR = randomGraphMyInteger(2, 3431, 0.2);
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexRandom4Less) {
  Graph<MyInteger> GR = randomGraphMyInteger(4, 3431, 0.2);
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(GR, 1);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 1);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexRandom8Less) {
  Graph<MyInteger> GR = randomGraphMyInteger(8, 3431, 0.2);
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexRandom16Less) {
  Graph<MyInteger> GR = randomGraphMyInteger(16, 3431, 0.2);
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(GR, 15);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 15);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexRandom32Less) {
  Graph<MyInteger> GR = randomGraphMyInteger(32, 3431, 0.2);
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexRandom64Less) {
  Graph<MyInteger> GR = randomGraphMyInteger(64, 3431, 0.2);
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexRandom128Less) {
  Graph<MyInteger> GR = randomGraphMyInteger(128, 3431, 0.2);
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

// Perform tests on same sized random graphs with more edges
TEST(DijkstraTest, distanceFrom0IndexRandom1More) {
  Graph<MyInteger> GR = randomGraphMyInteger(1, 3431, 0.8);
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexRandom2More) {
  Graph<MyInteger> GR = randomGraphMyInteger(2, 3431, 0.8);
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexRandom4More) {
  Graph<MyInteger> GR = randomGraphMyInteger(4, 3431, 0.8);
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexRandom8More) {
  Graph<MyInteger> GR = randomGraphMyInteger(8, 3431, 0.8);
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexRandom16More) {
  Graph<MyInteger> GR = randomGraphMyInteger(16, 3431, 0.8);
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexRandom32More) {
  Graph<MyInteger> GR = randomGraphMyInteger(32, 3431, 0.8);
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexRandom64More) {
  Graph<MyInteger> GR = randomGraphMyInteger(64, 3431, 0.8);
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0IndexRandom128More) {
  Graph<MyInteger> GR = randomGraphMyInteger(128, 3431, 0.8);
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceIndex<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

// singleSourceLazy tests /////////////////////////////////////////////////
// Test on int
TEST(DijkstraTest, distanceFrom0LazyTiny) {
  Graph<int> G {"tinyEWD.txt"};
  Graph<int> result =  singleSourceLazy<int>(G, 0);
  bool expected = correctShortestPathTree(G, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyMedium) {
  Graph<int> G {"mediumEWD.txt"};
  Graph<int> result =  singleSourceLazy<int>(G, 0);
  bool expected = correctShortestPathTree(G, result, 0);
  EXPECT_EQ(expected, true);
}

// Test on double
TEST(DijkstraTest, distanceFrom0LazyTinyDouble) {
  Graph<double> G {"tinyEWD.txt"};
  Graph<double> result =  singleSourceLazy<double>(G, 0);
  bool expected = correctShortestPathTree(G, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyMediumDouble) {
  Graph<double> G {"mediumEWD.txt"};
  Graph<double> result =  singleSourceLazy<double>(G, 0);
  bool expected = correctShortestPathTree(G, result, 0);
  EXPECT_EQ(expected, true);
}

// Test on MyInteger
TEST(DijkstraTest, distanceFrom0LazyTinyMyInt) {
  Graph<MyInteger> G {"tinyEWD.txt"};
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(G, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(G, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyMediumMyInt) {
  Graph<MyInteger> G {"mediumEWD.txt"};
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(G, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(G, result, 0);
  EXPECT_EQ(expected, true);
}

// Test from different sources
TEST(DijkstraTest, distanceFrom3LazyTinyMyInt) {
  Graph<MyInteger> G {"tinyEWD.txt"};
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(G, 3);
  bool expected = correctShortestPathTree(G, result, 3);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom3LazyMediumMyInt) {
  Graph<MyInteger> G {"mediumEWD.txt"};
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(G, 3);
  bool expected = correctShortestPathTree(G, result, 3);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom6LazyTinyMyInt) {
  Graph<MyInteger> G {"tinyEWD.txt"};
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(G, 6);
  bool expected = correctShortestPathTree(G, result, 6);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom8LazyMediumMyInt) {
  Graph<MyInteger> G {"mediumEWD.txt"};
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(G, 8);
  bool expected = correctShortestPathTree(G, result, 8);
  EXPECT_EQ(expected, true);
}

// Test on USA road
TEST(DijkstraTest, distanceFrom0LazyUSAMyInt) {
  Graph<MyInteger> G {"USA-road-d.NY.gr"};
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceLazy<MyInteger>(G, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(G, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom10LazyUSAMyInt) {
  Graph<MyInteger> G {"USA-road-d.NY.gr"};
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceLazy<MyInteger>(G, 10);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(G, result, 10);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom100LazyUSAMyInt) {
  Graph<MyInteger> G {"USA-road-d.NY.gr"};
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceLazy<MyInteger>(G, 100);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(G, result, 100);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom1000LazyUSAMyInt) {
  Graph<MyInteger> G {"USA-road-d.NY.gr"};
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceLazy<MyInteger>(G, 1000);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(G, result, 1000);
  EXPECT_EQ(expected, true);
}

// Test on smaller graphs with specific error prone characteristics
TEST(DijkstraTest, distanceFrom0LazyG1) {
  Graph<int> G1 {3};
  G1.addEdge(0, 1, 2);
  G1.addEdge(1, 0, 1);
  G1.addEdge(1, 2, 10);
  G1.addEdge(2, 1, 5);
  Graph<int> result =  singleSourceLazy<int>(G1, 0);
  bool expected = correctShortestPathTree(G1, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyG2) {
  Graph<double> G2 {3};
  G2.addEdge(0, 1, 2.1);
  G2.addEdge(1, 0, 1.4);
  G2.addEdge(1, 2, 10.0);
  G2.addEdge(2, 1, 5.89);
  Graph<double> result =  singleSourceLazy<double>(G2, 0);
  bool expected = correctShortestPathTree(G2, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyG3) {
  Graph<double> G3 {3};
  G3.addEdge(0, 1, 2.1);
  G3.addEdge(1, 0, 1.4);
  G3.addEdge(1, 2, 10.0);
  G3.addEdge(2, 1, 5.89);
  Graph<double> result =  singleSourceLazy<double>(G3, 0);
  bool expected = correctShortestPathTree(G3, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyG4) {
  Graph<double> G4 {2};
  G4.addEdge(0, 1, 2.1);
  G4.addEdge(1, 0, 1.4);
  Graph<double> result =  singleSourceLazy<double>(G4, 0);
  bool expected = correctShortestPathTree(G4, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyG5) {
  Graph<double> G5 {3};
  G5.addEdge(0, 1, 2.1);
  G5.addEdge(1, 0, 1.4);
  G5.addEdge(1, 2, 10.0);
  G5.addEdge(2, 1, 5.89);
  Graph<double> result =  singleSourceLazy<double>(G5, 0);
  bool expected = correctShortestPathTree(G5, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyG6) {
  Graph<double> G6 {4};
  G6.addEdge(0, 1, 2.1);
  G6.addEdge(1, 0, 1.4);
  G6.addEdge(1, 2, 10.0);
  G6.addEdge(2, 1, 5.89);
  G6.addEdge(3, 1, 1.27);
  Graph<double> result =  singleSourceLazy<double>(G6, 0);
  bool expected = correctShortestPathTree(G6, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyG7) {
  Graph<double> G7 {5};
  G7.addEdge(0, 2, 2.1);
  G7.addEdge(0, 4, 1.4);
  G7.addEdge(4, 1, 10.5);
  G7.addEdge(1, 3, 2.7);
  Graph<double> result =  singleSourceLazy<double>(G7, 0);
  bool expected = correctShortestPathTree(G7, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyG8) {
  Graph<double> G8 {6};
  G8.addEdge(0, 2, 2.1);
  G8.addEdge(0, 5, 1.4);
  G8.addEdge(5, 1, 10.5);
  G8.addEdge(1, 3, 2.7);
  Graph<double> result =  singleSourceLazy<double>(G8, 0);
  bool expected = correctShortestPathTree(G8, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyG9) {
  Graph<double> G9 {10};
  G9.addEdge(0, 2, 2.1);
  G9.addEdge(0, 4, 1.4);
  G9.addEdge(4, 1, 10.5);
  G9.addEdge(1, 3, 2.7);
  G9.addEdge(0, 9, 2.4);
  G9.addEdge(2, 7, 8.8);
  G9.addEdge(5, 6, 9.3);
  Graph<double> result =  singleSourceLazy<double>(G9, 0);
  bool expected = correctShortestPathTree(G9, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyG10) {
  Graph<double> G10 {7};
  G10.addEdge(0, 2, 2.1);
  G10.addEdge(0, 4, 1.4);
  G10.addEdge(1, 4, 10.5);
  G10.addEdge(1, 3, 2.7);
  G10.addEdge(0, 5, 2.4);
  G10.addEdge(2, 6, 8.8);
  Graph<double> result =  singleSourceLazy<double>(G10, 0);
  bool expected = correctShortestPathTree(G10, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyG11) {
  Graph<double> G11 {7};
  G11.addEdge(0, 2, 2.1);
  G11.addEdge(0, 4, 1.4);
  G11.addEdge(4, 1, 10.5);
  G11.addEdge(1, 3, 2.7);
  G11.addEdge(0, 5, 2.4);
  G11.addEdge(2, 6, 8.8);
  G11.addEdge(6, 5, 6.2);
  Graph<double> result =  singleSourceLazy<double>(G11, 0);
  bool expected = correctShortestPathTree(G11, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyG12) {
  Graph<int> G12 {4};
  G12.addEdge(3, 1, 4);
  G12.addEdge(3, 0, 1);
  G12.addEdge(3, 2, 2);
  G12.addEdge(1, 0, 9);
  Graph<int> result =  singleSourceLazy<int>(G12, 3);
  bool expected = correctShortestPathTree(G12, result, 3);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyG13) {
  Graph<int> G13 {7};
  G13.addEdge(0, 1, 2);
  G13.addEdge(0, 2, 1);
  G13.addEdge(1, 3, 1);
  G13.addEdge(1, 4, 5);
  G13.addEdge(2, 5, 3);
  G13.addEdge(2, 6, 2);
  Graph<int> result =  singleSourceLazy<int>(G13, 0);
  bool expected = correctShortestPathTree(G13, result, 0);
  EXPECT_EQ(expected, true);
}


TEST(DijkstraTest, distanceFrom0LazyG14) {
  Graph<int> G14 {7};
  G14.addEdge(0, 1, 1);
  G14.addEdge(0, 2, 1);
  G14.addEdge(1, 3, 1);
  G14.addEdge(1, 4, 1);
  G14.addEdge(2, 5, 1);
  G14.addEdge(2, 6, 1);
  Graph<int> result =  singleSourceLazy<int>(G14, 0);
  bool expected = correctShortestPathTree(G14, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom3LazyG15) {
  Graph<int> G15 {8};
  G15.addEdge(3, 2, 1);
  G15.addEdge(3, 4, 4);
  G15.addEdge(7, 5, 10);
  G15.addEdge(7, 6, 2);
  G15.addEdge(3, 7, 3);
  G15.addEdge(4, 1, 5);
  G15.addEdge(2, 0, 11);
  Graph<int> result =  singleSourceLazy<int>(G15, 3);
  bool expected = correctShortestPathTree(G15, result, 3);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyG16) {
  Graph<int> G16 {8};
  G16.addEdge(3, 2, 4);
  G16.addEdge(3, 4, 2);
  G16.addEdge(7, 5, 10);
  G16.addEdge(7, 6, 1);
  G16.addEdge(3, 7, 9);
  G16.addEdge(4, 1, 5);
  G16.addEdge(2, 0, 2);
  Graph<int> result =  singleSourceLazy<int>(G16, 0);
  bool expected = correctShortestPathTree(G16, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyG17) {
  Graph<int> G17(1);
  Graph<int> result =  singleSourceLazy<int>(G17, 0);
  bool expected = correctShortestPathTree(G17, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyG18) {
  Graph<int> G18(6);
  G18.addEdge(0, 1, 1);
  G18.addEdge(1, 2, 1);
  G18.addEdge(2, 3, 1);
  G18.addEdge(3, 4, 1);
  G18.addEdge(4, 5, 1);
  Graph<int> result =  singleSourceLazy<int>(G18, 0);
  bool expected = correctShortestPathTree(G18, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyG19) {
  Graph<int> G19(6);
  G19.addEdge(0, 1, 1);
  G19.addEdge(1, 2, 1);
  G19.addEdge(2, 3, 1);
  G19.addEdge(3, 4, 1);
  G19.addEdge(4, 5, 1);
  G19.addEdge(5, 0, 1);
  G19.addEdge(0, 3, 1);
  Graph<int> result =  singleSourceLazy<int>(G19, 0);
  bool expected = correctShortestPathTree(G19, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyG20) {
  Graph<int> G20(4);
  G20.addEdge(0, 1, 4); 
  G20.addEdge(0, 2, 5); 
  G20.addEdge(2, 1, 3); 
  G20.addEdge(1, 3, 1); 
  Graph<int> result =  singleSourceLazy<int>(G20, 0);
  bool expected = correctShortestPathTree(G20, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyG21) {
  Graph<int> G21(6);
  G21.addEdge(0, 1, 1); 
  G21.addEdge(0, 2, 1); 
  G21.addEdge(0, 3, 1); 
  G21.addEdge(0, 4, 1); 
  G21.addEdge(0, 5, 1); 
  G21.addEdge(2, 4, 1); 
  G21.addEdge(4, 5, 1); 
  Graph<int> result =  singleSourceLazy<int>(G21, 0);
  bool expected = correctShortestPathTree(G21, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyG22) {
  Graph<int> G22(5);
  G22.addEdge(0, 1, 4);
  G22.addEdge(1, 2, 2);
  G22.addEdge(2, 3, 2);
  G22.addEdge(3, 4, 1);
  G22.addEdge(0, 3, 2);
  Graph<int> result =  singleSourceLazy<int>(G22, 0);
  bool expected = correctShortestPathTree(G22, result, 0);
  EXPECT_EQ(expected, true);
}

// Test on random graphs of varying sizes
TEST(DijkstraTest, distanceFrom0LazyRandom1) {
  Graph<MyInteger> GR = randomGraphMyInteger(1, 3431);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyRandom2) {
  Graph<MyInteger> GR = randomGraphMyInteger(2, 3431);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyRandom4) {
  Graph<MyInteger> GR = randomGraphMyInteger(4, 3431);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyRandom8) {
  Graph<MyInteger> GR = randomGraphMyInteger(8, 3431);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyRandom16) {
  Graph<MyInteger> GR = randomGraphMyInteger(16, 3431);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyRandom32) {
  Graph<MyInteger> GR = randomGraphMyInteger(32, 3431);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyRandom64) {
  Graph<MyInteger> GR = randomGraphMyInteger(64, 3431);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyRandom128) {
  Graph<MyInteger> GR = randomGraphMyInteger(128, 3431);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

// Perform tests on same sized random graphs with less edges
TEST(DijkstraTest, distanceFrom0LazyRandom1Less) {
  Graph<MyInteger> GR = randomGraphMyInteger(1, 3431, 0.2);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyRandom2Less) {
  Graph<MyInteger> GR = randomGraphMyInteger(2, 3431, 0.2);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyRandom4Less) {
  Graph<MyInteger> GR = randomGraphMyInteger(4, 3431, 0.2);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(GR, 1);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 1);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyRandom8Less) {
  Graph<MyInteger> GR = randomGraphMyInteger(8, 3431, 0.2);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyRandom16Less) {
  Graph<MyInteger> GR = randomGraphMyInteger(16, 3431, 0.2);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(GR, 15);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 15);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyRandom32Less) {
  Graph<MyInteger> GR = randomGraphMyInteger(32, 3431, 0.2);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyRandom64Less) {
  Graph<MyInteger> GR = randomGraphMyInteger(64, 3431, 0.2);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyRandom128Less) {
  Graph<MyInteger> GR = randomGraphMyInteger(128, 3431, 0.2);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

// Perform tests on same sized random graphs with more edges
TEST(DijkstraTest, distanceFrom0LazyRandom1More) {
  Graph<MyInteger> GR = randomGraphMyInteger(1, 3431, 0.8);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyRandom2More) {
  Graph<MyInteger> GR = randomGraphMyInteger(2, 3431, 0.8);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyRandom4More) {
  Graph<MyInteger> GR = randomGraphMyInteger(4, 3431, 0.8);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyRandom8More) {
  Graph<MyInteger> GR = randomGraphMyInteger(8, 3431, 0.8);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyRandom16More) {
  Graph<MyInteger> GR = randomGraphMyInteger(16, 3431, 0.8);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyRandom32More) {
  Graph<MyInteger> GR = randomGraphMyInteger(32, 3431, 0.8);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyRandom64More) {
  Graph<MyInteger> GR = randomGraphMyInteger(64, 3431, 0.8);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0LazyRandom128More) {
  Graph<MyInteger> GR = randomGraphMyInteger(128, 3431, 0.8);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceLazy<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

// singleSourceSet tests ////////////////////////////////////////////////
// Test on int
TEST(DijkstraTest, distanceFrom0SetTiny) {
  Graph<int> G {"tinyEWD.txt"};
  Graph<int> result =  singleSourceSet<int>(G, 0);
  bool expected = correctShortestPathTree(G, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetMedium) {
  Graph<int> G {"mediumEWD.txt"};
  Graph<int> result =  singleSourceSet<int>(G, 0);
  bool expected = correctShortestPathTree(G, result, 0);
  EXPECT_EQ(expected, true);
}

// Test on double
TEST(DijkstraTest, distanceFrom0SetTinyDouble) {
  Graph<double> G {"tinyEWD.txt"};
  Graph<double> result =  singleSourceSet<double>(G, 0);
  bool expected = correctShortestPathTree(G, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetMediumDouble) {
  Graph<double> G {"mediumEWD.txt"};
  Graph<double> result =  singleSourceSet<double>(G, 0);
  bool expected = correctShortestPathTree(G, result, 0);
  EXPECT_EQ(expected, true);
}

// Test on MyInteger
TEST(DijkstraTest, distanceFrom0SetTinyMyInt) {
  Graph<MyInteger> G {"tinyEWD.txt"};
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(G, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(G, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetMediumMyInt) {
  Graph<MyInteger> G {"mediumEWD.txt"};
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(G, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(G, result, 0);
  EXPECT_EQ(expected, true);
}

// Test from different sources
TEST(DijkstraTest, distanceFrom3SetTinyMyInt) {
  Graph<MyInteger> G {"tinyEWD.txt"};
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(G, 3);
  bool expected = correctShortestPathTree(G, result, 3);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom3SetMediumMyInt) {
  Graph<MyInteger> G {"mediumEWD.txt"};
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(G, 3);
  bool expected = correctShortestPathTree(G, result, 3);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom6SetTinyMyInt) {
  Graph<MyInteger> G {"tinyEWD.txt"};
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(G, 6);
  bool expected = correctShortestPathTree(G, result, 6);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom8SetMediumMyInt) {
  Graph<MyInteger> G {"mediumEWD.txt"};
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(G, 8);
  bool expected = correctShortestPathTree(G, result, 8);
  EXPECT_EQ(expected, true);
}

// Test on USA road
TEST(DijkstraTest, distanceFrom0SetUSAMyInt) {
  Graph<MyInteger> G {"USA-road-d.NY.gr"};
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceSet<MyInteger>(G, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(G, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom10SetUSAMyInt) {
  Graph<MyInteger> G {"USA-road-d.NY.gr"};
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceSet<MyInteger>(G, 10);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(G, result, 10);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom100SetUSAMyInt) {
  Graph<MyInteger> G {"USA-road-d.NY.gr"};
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceSet<MyInteger>(G, 100);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(G, result, 100);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom1000SetUSAMyInt) {
  Graph<MyInteger> G {"USA-road-d.NY.gr"};
  MyInteger::clearCounts();
  Graph<MyInteger> result = singleSourceSet<MyInteger>(G, 1000);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(G, result, 1000);
  EXPECT_EQ(expected, true);
}

// Test on smaller graphs with specific error prone characteristics
TEST(DijkstraTest, distanceFrom0SetG1) {
  Graph<int> G1 {3};
  G1.addEdge(0, 1, 2);
  G1.addEdge(1, 0, 1);
  G1.addEdge(1, 2, 10);
  G1.addEdge(2, 1, 5);
  Graph<int> result =  singleSourceSet<int>(G1, 0);
  bool expected = correctShortestPathTree(G1, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetG2) {
  Graph<double> G2 {3};
  G2.addEdge(0, 1, 2.1);
  G2.addEdge(1, 0, 1.4);
  G2.addEdge(1, 2, 10.0);
  G2.addEdge(2, 1, 5.89);
  Graph<double> result =  singleSourceSet<double>(G2, 0);
  bool expected = correctShortestPathTree(G2, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetG3) {
  Graph<double> G3 {3};
  G3.addEdge(0, 1, 2.1);
  G3.addEdge(1, 0, 1.4);
  G3.addEdge(1, 2, 10.0);
  G3.addEdge(2, 1, 5.89);
  Graph<double> result =  singleSourceSet<double>(G3, 0);
  bool expected = correctShortestPathTree(G3, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetG4) {
  Graph<double> G4 {2};
  G4.addEdge(0, 1, 2.1);
  G4.addEdge(1, 0, 1.4);
  Graph<double> result =  singleSourceSet<double>(G4, 0);
  bool expected = correctShortestPathTree(G4, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetG5) {
  Graph<double> G5 {3};
  G5.addEdge(0, 1, 2.1);
  G5.addEdge(1, 0, 1.4);
  G5.addEdge(1, 2, 10.0);
  G5.addEdge(2, 1, 5.89);
  Graph<double> result =  singleSourceSet<double>(G5, 0);
  bool expected = correctShortestPathTree(G5, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetG6) {
  Graph<double> G6 {4};
  G6.addEdge(0, 1, 2.1);
  G6.addEdge(1, 0, 1.4);
  G6.addEdge(1, 2, 10.0);
  G6.addEdge(2, 1, 5.89);
  G6.addEdge(3, 1, 1.27);
  Graph<double> result =  singleSourceSet<double>(G6, 0);
  bool expected = correctShortestPathTree(G6, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetG7) {
  Graph<double> G7 {5};
  G7.addEdge(0, 2, 2.1);
  G7.addEdge(0, 4, 1.4);
  G7.addEdge(4, 1, 10.5);
  G7.addEdge(1, 3, 2.7);
  Graph<double> result =  singleSourceSet<double>(G7, 0);
  bool expected = correctShortestPathTree(G7, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetG8) {
  Graph<double> G8 {6};
  G8.addEdge(0, 2, 2.1);
  G8.addEdge(0, 5, 1.4);
  G8.addEdge(5, 1, 10.5);
  G8.addEdge(1, 3, 2.7);
  Graph<double> result =  singleSourceSet<double>(G8, 0);
  bool expected = correctShortestPathTree(G8, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetG9) {
  Graph<double> G9 {10};
  G9.addEdge(0, 2, 2.1);
  G9.addEdge(0, 4, 1.4);
  G9.addEdge(4, 1, 10.5);
  G9.addEdge(1, 3, 2.7);
  G9.addEdge(0, 9, 2.4);
  G9.addEdge(2, 7, 8.8);
  G9.addEdge(5, 6, 9.3);
  Graph<double> result =  singleSourceSet<double>(G9, 0);
  bool expected = correctShortestPathTree(G9, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetG10) {
  Graph<double> G10 {7};
  G10.addEdge(0, 2, 2.1);
  G10.addEdge(0, 4, 1.4);
  G10.addEdge(1, 4, 10.5);
  G10.addEdge(1, 3, 2.7);
  G10.addEdge(0, 5, 2.4);
  G10.addEdge(2, 6, 8.8);
  Graph<double> result =  singleSourceSet<double>(G10, 0);
  bool expected = correctShortestPathTree(G10, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetG11) {
  Graph<double> G11 {7};
  G11.addEdge(0, 2, 2.1);
  G11.addEdge(0, 4, 1.4);
  G11.addEdge(4, 1, 10.5);
  G11.addEdge(1, 3, 2.7);
  G11.addEdge(0, 5, 2.4);
  G11.addEdge(2, 6, 8.8);
  G11.addEdge(6, 5, 6.2);
  Graph<double> result =  singleSourceSet<double>(G11, 0);
  bool expected = correctShortestPathTree(G11, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetG12) {
  Graph<int> G12 {4};
  G12.addEdge(3, 1, 4);
  G12.addEdge(3, 0, 1);
  G12.addEdge(3, 2, 2);
  G12.addEdge(1, 0, 9);
  Graph<int> result =  singleSourceSet<int>(G12, 3);
  bool expected = correctShortestPathTree(G12, result, 3);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetG13) {
  Graph<int> G13 {7};
  G13.addEdge(0, 1, 2);
  G13.addEdge(0, 2, 1);
  G13.addEdge(1, 3, 1);
  G13.addEdge(1, 4, 5);
  G13.addEdge(2, 5, 3);
  G13.addEdge(2, 6, 2);
  Graph<int> result =  singleSourceSet<int>(G13, 0);
  bool expected = correctShortestPathTree(G13, result, 0);
  EXPECT_EQ(expected, true);
}


TEST(DijkstraTest, distanceFrom0SetG14) {
  Graph<int> G14 {7};
  G14.addEdge(0, 1, 1);
  G14.addEdge(0, 2, 1);
  G14.addEdge(1, 3, 1);
  G14.addEdge(1, 4, 1);
  G14.addEdge(2, 5, 1);
  G14.addEdge(2, 6, 1);
  Graph<int> result =  singleSourceSet<int>(G14, 0);
  bool expected = correctShortestPathTree(G14, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom3SetG15) {
  Graph<int> G15 {8};
  G15.addEdge(3, 2, 1);
  G15.addEdge(3, 4, 4);
  G15.addEdge(7, 5, 10);
  G15.addEdge(7, 6, 2);
  G15.addEdge(3, 7, 3);
  G15.addEdge(4, 1, 5);
  G15.addEdge(2, 0, 11);
  Graph<int> result =  singleSourceSet<int>(G15, 3);
  bool expected = correctShortestPathTree(G15, result, 3);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetG16) {
  Graph<int> G16 {8};
  G16.addEdge(3, 2, 4);
  G16.addEdge(3, 4, 2);
  G16.addEdge(7, 5, 10);
  G16.addEdge(7, 6, 1);
  G16.addEdge(3, 7, 9);
  G16.addEdge(4, 1, 5);
  G16.addEdge(2, 0, 2);
  Graph<int> result =  singleSourceSet<int>(G16, 0);
  bool expected = correctShortestPathTree(G16, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetG17) {
  Graph<int> G17(1);
  Graph<int> result =  singleSourceSet<int>(G17, 0);
  bool expected = correctShortestPathTree(G17, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetG18) {
  Graph<int> G18(6);
  G18.addEdge(0, 1, 1);
  G18.addEdge(1, 2, 1);
  G18.addEdge(2, 3, 1);
  G18.addEdge(3, 4, 1);
  G18.addEdge(4, 5, 1);
  Graph<int> result =  singleSourceSet<int>(G18, 0);
  bool expected = correctShortestPathTree(G18, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetG19) {
  Graph<int> G19(6);
  G19.addEdge(0, 1, 1);
  G19.addEdge(1, 2, 1);
  G19.addEdge(2, 3, 1);
  G19.addEdge(3, 4, 1);
  G19.addEdge(4, 5, 1);
  G19.addEdge(5, 0, 1);
  G19.addEdge(0, 3, 1);
  Graph<int> result =  singleSourceSet<int>(G19, 0);
  bool expected = correctShortestPathTree(G19, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetG20) {
  Graph<int> G20(4);
  G20.addEdge(0, 1, 4); 
  G20.addEdge(0, 2, 5); 
  G20.addEdge(2, 1, 3); 
  G20.addEdge(1, 3, 1); 
  Graph<int> result =  singleSourceSet<int>(G20, 0);
  bool expected = correctShortestPathTree(G20, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetG21) {
  Graph<int> G21(6);
  G21.addEdge(0, 1, 1); 
  G21.addEdge(0, 2, 1); 
  G21.addEdge(0, 3, 1); 
  G21.addEdge(0, 4, 1); 
  G21.addEdge(0, 5, 1); 
  G21.addEdge(2, 4, 1); 
  G21.addEdge(4, 5, 1); 
  Graph<int> result =  singleSourceSet<int>(G21, 0);
  bool expected = correctShortestPathTree(G21, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetG22) {
  Graph<int> G22(5);
  G22.addEdge(0, 1, 4);
  G22.addEdge(1, 2, 2);
  G22.addEdge(2, 3, 2);
  G22.addEdge(3, 4, 1);
  G22.addEdge(0, 3, 2);
  Graph<int> result =  singleSourceSet<int>(G22, 0);
  bool expected = correctShortestPathTree(G22, result, 0);
  EXPECT_EQ(expected, true);
}

// Test on random graphs of varying sizes
TEST(DijkstraTest, distanceFrom0SetRandom1) {
  Graph<MyInteger> GR = randomGraphMyInteger(1, 3431);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetRandom2) {
  Graph<MyInteger> GR = randomGraphMyInteger(2, 3431);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetRandom4) {
  Graph<MyInteger> GR = randomGraphMyInteger(4, 3431);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetRandom8) {
  Graph<MyInteger> GR = randomGraphMyInteger(8, 3431);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetRandom16) {
  Graph<MyInteger> GR = randomGraphMyInteger(16, 3431);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetRandom32) {
  Graph<MyInteger> GR = randomGraphMyInteger(32, 3431);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetRandom64) {
  Graph<MyInteger> GR = randomGraphMyInteger(64, 3431);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetRandom128) {
  Graph<MyInteger> GR = randomGraphMyInteger(128, 3431);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

// Perform tests on same sized random graphs with less edges
TEST(DijkstraTest, distanceFrom0SetRandom1Less) {
  Graph<MyInteger> GR = randomGraphMyInteger(1, 3431, 0.2);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetRandom2Less) {
  Graph<MyInteger> GR = randomGraphMyInteger(2, 3431, 0.2);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetRandom4Less) {
  Graph<MyInteger> GR = randomGraphMyInteger(4, 3431, 0.2);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(GR, 1);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 1);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetRandom8Less) {
  Graph<MyInteger> GR = randomGraphMyInteger(8, 3431, 0.2);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetRandom16Less) {
  Graph<MyInteger> GR = randomGraphMyInteger(16, 3431, 0.2);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(GR, 15);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 15);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetRandom32Less) {
  Graph<MyInteger> GR = randomGraphMyInteger(32, 3431, 0.2);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetRandom64Less) {
  Graph<MyInteger> GR = randomGraphMyInteger(64, 3431, 0.2);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetRandom128Less) {
  Graph<MyInteger> GR = randomGraphMyInteger(128, 3431, 0.2);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

// Perform tests on same sized random graphs with more edges
TEST(DijkstraTest, distanceFrom0SetRandom1More) {
  Graph<MyInteger> GR = randomGraphMyInteger(1, 3431, 0.8);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetRandom2More) {
  Graph<MyInteger> GR = randomGraphMyInteger(2, 3431, 0.8);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetRandom4More) {
  Graph<MyInteger> GR = randomGraphMyInteger(4, 3431, 0.8);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetRandom8More) {
  Graph<MyInteger> GR = randomGraphMyInteger(8, 3431, 0.8);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetRandom16More) {
  Graph<MyInteger> GR = randomGraphMyInteger(16, 3431, 0.8);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetRandom32More) {
  Graph<MyInteger> GR = randomGraphMyInteger(32, 3431, 0.8);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetRandom64More) {
  Graph<MyInteger> GR = randomGraphMyInteger(64, 3431, 0.8);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

TEST(DijkstraTest, distanceFrom0SetRandom128More) {
  Graph<MyInteger> GR = randomGraphMyInteger(128, 3431, 0.8);
  MyInteger::clearCounts();
  Graph<MyInteger> result =  singleSourceSet<MyInteger>(GR, 0);
  MyInteger::printCounts();
  bool expected = correctShortestPathTree(GR, result, 0);
  EXPECT_EQ(expected, true);
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
