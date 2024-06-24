#ifndef GRAPH_HPP_
#define GRAPH_HPP_

#include <iostream>
#include <fstream>
#include <utility>
#include <functional>
#include <vector>
#include <string>
#include <queue>
#include <set>
#include <unordered_map>
#include <limits>
#include "my_integer.hpp"
#include <stack>
#include "main.cpp"

// Sources
// https://bradfieldcs.com/algos/graphs/dijkstras-algorithm/#:~:text=Dijkstra's%20algorithm%20uses%20a%20priority,of%20vertices%20sorted%20by%20distance.
// https://stackoverflow.com/questions/12782431/relaxation-of-an-edge-in-dijkstras-algorithm
// https://www.javatpoint.com/cpp-dijkstra-algorithm-using-priority-queue

template <typename T>
class Graph {
 private:
  std::vector<std::unordered_map<int, T> > adjList {};
  int numVertices {};

 public:
  // empty graph with N vertices
  explicit Graph(int N);

  // construct graph from edge list in filename
  explicit Graph(const std::string& filename);

  // add an edge directed from vertex i to vertex j with given weight
  void addEdge(int i, int j, T weight);

  // removes edge from vertex i to vertex j
  void removeEdge(int i, int j);

  // is there an edge from vertex i to vertex j?
  bool isEdge(int i, int j) const;

  // return weight of edge from i to j
  // will throw an exception if there is no edge from i to j
  T getEdgeWeight(int i, int j) const;

  // returns number of vertices in the graph
  int size() const;

  // alias a const iterator to our adjacency list type to iterator
  using iterator = 
  typename std::vector<std::unordered_map<int, T> >::const_iterator;

  // cbegin returns const iterator pointing to first element of adjList
  iterator begin() const {
    return adjList.cbegin();
  }

  iterator end() const {
    return adjList.cend();
  }

  // return iterator to a particular vertex
  iterator neighbours(int a) const {
    return adjList.begin() + a;
  }
};

template <typename T>
Graph<T>::Graph(int N) : adjList(N), numVertices {N} {}

template <typename T>
Graph<T>::Graph(const std::string& inputFile) {
  std::ifstream infile {inputFile};
  if (!infile) {
    std::cerr << inputFile << " could not be opened\n";
    return;
  }
  // first line has number of vertices
  infile >> numVertices;
  adjList.resize(numVertices);
  int i {};
  int j {};
  double weight {};
  // assume each remaining line is of form
  // origin dest weight
  while (infile >> i >> j >> weight) {
    addEdge(i, j, static_cast<T>(weight));
  }
}

template <typename T>
int Graph<T>::size() const {
  return numVertices;
}

template <typename T>
void Graph<T>::addEdge(int i, int j, T weight) {
  if (i < 0 or i >= numVertices or j < 0 or j >= numVertices) {
    throw std::out_of_range("invalid vertex number");
  }
  adjList[i].insert({j, weight});
}

template <typename T>
void Graph<T>::removeEdge(int i, int j) {
  // check if i and j are valid
  if (i >= 0 && i < numVertices && j >= 0 && j < numVertices) {
    adjList[i].erase(j);
  }
}

template <typename T>
bool Graph<T>::isEdge(int i, int j) const {
  if (i >= 0 && i < numVertices && j >= 0 && j < numVertices) {
    return adjList.at(i).contains(j);
  }
  return false;
}

template <typename T>
T Graph<T>::getEdgeWeight(int i, int j) const {
  return adjList.at(i).at(j);
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const Graph<T>& G) {
  for (int i = 0; i < G.size(); ++i) {
    out << i << ':';
    for (const auto& edge : *(G.neighbours(i))) {
      out << " (" << i << ", " << edge.first << ")[" << edge.second << ']';
    }
    out << '\n';
  }
  return out;
}

// End of functions from Graph class

template <typename T>
T infinity() {
  if (std::numeric_limits<T>::has_infinity) {
    return std::numeric_limits<T>::infinity();
  } else {
    return std::numeric_limits<T>::max();
  }
}

template <typename T>
class IndexPriorityQueue {
 private:
  std::vector<T> priorities;
  std::vector<int> priorityQueue;
  std::vector<int> indexToPosition;
  
  int size_ = 0;
  int maxSize_;

 public:
  explicit IndexPriorityQueue(int N);
  void push(const T& priority, int index);
  void pop();
  void erase(int index);
  bool contains(int index) const;
  void changeKey(const T& priority, int index);
  std::pair<T, int> top() const;
  bool empty() const;
  int size() const;

 private:
  void swim(int i);
  void sink(int i);
  void swap(int i, int j);
};

// Useful helper functions
int leftChild(int i) {
  return 2 * i;
}

int rightChild(int i) {
  return 2 * i + 1;
}

int parent(int i) {
  return i / 2;
}

// IndexPriorityQueue member functions
// Initialize with maximum size N
template <typename T>
IndexPriorityQueue<T>::IndexPriorityQueue(int N) : maxSize_(N) {
  priorities.resize(N);
  priorityQueue.resize(N + 1); // 1-based index for the heap
  indexToPosition.resize(N, -1);
}

template <typename T>
bool IndexPriorityQueue<T>::empty() const {
  return size_ == 0;
}

template <typename T>
int IndexPriorityQueue<T>::size() const {
  return size_;
}

// Push a new element with priority and index
template <typename T>
void IndexPriorityQueue<T>::push(const T& priority, int index) {
  if (contains(index)) {
    return;
  }
  priorities.at(index) = priority;
  priorityQueue.at(++size_) = index;
  indexToPosition.at(index) = size_;
  swim(size_);
}

// Pop the element with the lowest priority
template <typename T>
void IndexPriorityQueue<T>::pop() {
  if (empty()) {
    throw std::out_of_range("Priority queue is empty");
  }
  int minIndex = priorityQueue.at(1);
  swap(1, size_--);
  sink(1);
  indexToPosition.at(minIndex) = -1;
}

// Erase the element with the given index
template <typename T>
void IndexPriorityQueue<T>::erase(int index) {
  if (!contains(index)) {
    return;
  }
  int pos = indexToPosition.at(index);
  swap(pos, size_--);
  swim(pos);
  sink(pos);
  indexToPosition.at(index) = -1;
}

// Return the element with the minimum priority and its index
template <typename T>
std::pair<T, int> IndexPriorityQueue<T>::top() const {
  if (empty()) {
    throw std::out_of_range("Priority queue is empty");
  }
  int minIndex = priorityQueue.at(1);
  return {priorities.at(minIndex), minIndex};
}

// Change the priority of the element with the given index
template <typename T>
void IndexPriorityQueue<T>::changeKey(const T& priority, int index) {
  if (!contains(index)) {
    push(priority, index);
  } else {
    priorities.at(index) = priority;
    int pos = indexToPosition.at(index);
    swim(pos);
    sink(pos);
  }
}

// Check if the element with the given index is in the priority queue
template <typename T>
bool IndexPriorityQueue<T>::contains(int index) const {
  if (index < 0 || index >= maxSize_) {
    return false;
  }
  return indexToPosition.at(index) != -1;
}

// Rearrange order by swimming up
template <typename T>
void IndexPriorityQueue<T>::swim(int i) {
  while (i > 1 && priorities.at(priorityQueue.at(i)) < priorities.at(priorityQueue.at(parent(i)))) {
    swap(i, parent(i));
    i = parent(i);
  }
}

// Maintain order by sinking down
template <typename T>
void IndexPriorityQueue<T>::sink(int i) {
  while (leftChild(i) <= size_) {
    int j = leftChild(i);
    if (j < size_ && priorities.at(priorityQueue.at(j)) > priorities.at(priorityQueue.at(j + 1))) {
      j++;
    }
    if (priorities.at(priorityQueue.at(i)) <= priorities.at(priorityQueue.at(j))) {
      break;
    }
    swap(i, j);
    i = j;
  }
}

template <typename T>
void IndexPriorityQueue<T>::swap(int i, int j) {
  std::swap(priorityQueue.at(i), priorityQueue.at(j));
  indexToPosition.at(priorityQueue.at(i)) = i;
  indexToPosition.at(priorityQueue.at(j)) = j;
}

// Implement your solution using an index priority queue here
template <typename T>
Graph<T> singleSourceIndex(const Graph<T>& G, int source) {
  int n = G.size(); // get the number of vertices
  std::vector<T> bestDistanceTo(n, infinity<T>()); // store best distance to each vertex from source
  std::vector<int> parent(n, -1); // store the parent of each vertex
  std::vector<T> edgeWeight(n, T{}); // store the weights of the edges in the shortest path tree
  bestDistanceTo[source] = T{}; 
  IndexPriorityQueue<T> pq(n); // initialize pq with number of vertic
  pq.push(T{}, source);

  while (!pq.empty()) {
    auto [dist, current] = pq.top();
    pq.pop();

    // relax all edges from the current vertex
    for (const auto& [neighbour, weight] : *(G.neighbours(current))) {
      // calculate distance to neighbour via current vertex
      T distanceViaCurrent = bestDistanceTo[current] + weight;

      // if shorter path found, update best distance to the neighbour
      if (bestDistanceTo[neighbour] > distanceViaCurrent) {
        bestDistanceTo[neighbour] = distanceViaCurrent;
        parent[neighbour] = current; // update parent of the neighbour to current
        edgeWeight[neighbour] = weight; // store the weight of this edge
        if (pq.contains(neighbour)) { // if neighbour is in pq, update its priority
          pq.changeKey(distanceViaCurrent, neighbour);
        } else {
          pq.push(distanceViaCurrent, neighbour); // or push neighbour into pq with new distance if not found
        }
      }
    }
  }
  // shortest path tree
  Graph<T> shortestPaths(n);

  // add edges to shortest path tree based on parent relationship
  for (int v = 0; v < n; ++v) {
    if (parent[v] != -1) {
      shortestPaths.addEdge(parent[v], v, edgeWeight[v]);
    }
  }
  return shortestPaths;
}

// Implement your lazy solution using std::priority_queue here
template <typename T>
Graph<T> singleSourceLazy(const Graph<T>& G, int source) {
  using distanceAndVertex = std::pair<T, int>;
  using minPQ = std::priority_queue<distanceAndVertex, std::vector<distanceAndVertex>, std::greater<distanceAndVertex>>; // start with the smallest 
  int n = G.size();

  std::vector<T> bestDistanceTo(n, infinity<T>()); // vector to store best known distance to each vertex
  std::vector<int> parent (n, -1); 
  std::vector<bool> visited (n, false);
  std::vector<T> edgeWeight(n, T{}); // vector to store the weights of the edges

  minPQ pq;
  bestDistanceTo[source] = T{}; //dist to source is 0
  pq.push({T{}, source}); // push source into the pq

  while(!pq.empty()){
    auto [dist, current] = pq.top();
    pq.pop();
    if(visited[current]){
      continue;
    }
    visited[current] = true;

    for(const auto& [neighbour, weight] : *(G.neighbours(current))){ // iterate over neighbours of the current vertex
      T distanceViaCurrent = bestDistanceTo[current] + weight; // calculate distance to the neighbour via the current vertex
      if(bestDistanceTo[neighbour] > distanceViaCurrent){ //if shorter path found update best known distance to the neighbour
        bestDistanceTo[neighbour] = distanceViaCurrent;
        parent[neighbour] = current; //parent of the neighbour is the current vertex in this case
        edgeWeight[neighbour] = weight; // store the weight of this edge
        pq.push({distanceViaCurrent, neighbour}); //push neighbour into the priority queue with updated distance
      }
    }
  }
  Graph<T> shortestPaths(n); //store shortest path tree
  for(int v = 0; v < n; ++v){
    if (parent[v] != -1){
      shortestPaths.addEdge(parent[v], v, edgeWeight[v]); // add edges based on the established parent relationships
    }
  }
  return shortestPaths;
}

// put your "best" solution here
// this is the one we will use for performance testing
template <typename T>
Graph<T> singleSourceShortestPaths(const Graph<T>& G, int source) {
  int n = G.size(); // get the number of vertices
  std::vector<T> bestDistanceTo(n, infinity<T>()); // store best distance to each vertex from source
  std::vector<int> parent(n, -1); // store the parent of each vertex
  std::vector<T> edgeWeight(n, T{}); // store the weights of the edges
  bestDistanceTo[source] = T{}; 
  IndexPriorityQueue<T> pq(n); // initialize pq with number of vertices
  pq.push(T{}, source);

  while (!pq.empty()) {
    auto [dist, current] = pq.top();
    pq.pop();

    // relax outgoing edges from the current vertex
    for (const auto& [neighbour, weight] : *(G.neighbours(current))) {
      // calculate distance to neighbour via current vertex
      T distanceViaCurrent = bestDistanceTo[current] + weight;

      // if shorter path found, update best distance to the neighbour
      if (bestDistanceTo[neighbour] > distanceViaCurrent) {
        bestDistanceTo[neighbour] = distanceViaCurrent;
        parent[neighbour] = current; // update parent of the neighbour to current
        edgeWeight[neighbour] = weight; // store the weight of this edge
        if (pq.contains(neighbour)) { // if neighbour is in pq, update its priority
          pq.changeKey(distanceViaCurrent, neighbour);
        } else {
          pq.push(distanceViaCurrent, neighbour); // or push neighbour into pq with new distance
        }
      }
    }
  }
  // shortest path tree
  Graph<T> shortestPaths(n);

  // add edges to shortest path tree based on parent relationship
  for (int v = 0; v < n; ++v) {
    if (parent[v] != -1) {
      shortestPaths.addEdge(parent[v], v, edgeWeight[v]);
    }
  }
  return shortestPaths;
}

#endif      // GRAPH_HPP_
