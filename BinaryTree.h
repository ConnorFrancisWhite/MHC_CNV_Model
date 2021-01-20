#ifndef _BINARYTREE
#define _BINARYTREE
using namespace std;

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "Datum.h"

class BinaryTree
{
  Datum data;

public:
  BinaryTree *L;
  BinaryTree *R;
  BinaryTree();
  BinaryTree(Datum d);
  ~BinaryTree();
  void build(BinaryTree *T, vector<Datum> a);
  void print(BinaryTree *T, int level, ofstream &out);
  Datum search(BinaryTree *T, double x, int level);
};

#endif

