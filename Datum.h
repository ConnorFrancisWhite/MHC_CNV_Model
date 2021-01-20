#ifndef _DATUM
#define _DATUM

// Define a simple data container class which stores an int and a double
class Datum
{
public:
  int n;
  double x;
  Datum();
  // Define a constructor which assigns values to n and x when a data object is 
  // instantiated
  Datum(int N, double X);
};

#endif

