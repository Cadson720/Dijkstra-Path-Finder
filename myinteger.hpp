#ifndef MY_INTEGER_HPP_
#define MY_INTEGER_HPP_

#include <utility>

struct MyInteger {
  int value {};
  inline static int constructorCount = 0;
  inline static int copyCount = 0;
  inline static int assignmentCount = 0;
  inline static int equalityCount = 0;
  inline static int comparisonCount = 0;

  explicit MyInteger(int initialValue = 0) : value {initialValue} {
    ++constructorCount;
  }
  
  ~MyInteger() {}

  MyInteger(const MyInteger& other) {
    ++copyCount;
    value = other.value;
  }

  MyInteger& operator=(const MyInteger& other) {
    ++assignmentCount;
    value = other.value;
    return *this;
  }

  static void printCounts() {
    std::cout << "constructor count is " << constructorCount << '\n';
    std::cout << "copy count is " << copyCount << '\n';
    std::cout << "assignment count is " << assignmentCount << '\n';
    std::cout << "comparison count is " << comparisonCount << '\n';
  }

  static void clearCounts() {
    constructorCount = 0;
    copyCount = 0;
    assignmentCount = 0;
    equalityCount = 0;
    comparisonCount = 0;
  }

  friend bool operator==(const MyInteger& x, const MyInteger& y) {
    ++equalityCount;
    return x.value == y.value;
  }

  friend bool operator!=(const MyInteger& x, const MyInteger& y) {
    return not (x == y);
  }

  friend auto operator<=>(const MyInteger& x, const MyInteger& y) {
    ++comparisonCount;
    return x.value <=> y.value;
  }

  friend std::ostream& operator<<(std::ostream& out, const MyInteger& a) {
    out << a.value;
    return out;
  }
};

MyInteger operator+(const MyInteger& x, const MyInteger& y) {
  return MyInteger {x.value + y.value};
}

#endif    // MY_INTEGER_HPP_
