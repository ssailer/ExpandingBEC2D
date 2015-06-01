#include <inttypes.h>
#include <iostream>
#include <stdio.h>
#include <python2.7/Python.h>

using namespace std;

int main(int argc, char *argv[])
{	
	setenv("PYTHONPATH",".",1);
    // Initialize the Python Interpreter
    Py_Initialize();

  	PyRun_SimpleString("import hello");
  	PyRun_SimpleString("hello.hello()");

    // Finish the Python Interpreter
    Py_Finalize();
    return 0;
}