#define PY_SSIZE_T_CLEAN
#include <Python.h>

static PyObject * access_system(PyObject *self, PyObject *args) 
{
   const char *info;
   double value;
   double *p;

   if (!PyArg_ParseTuple(args, "s", &info))
      return NULL;
   if (!PyArg_ParseTuple(args, "d", &value))
      return NULL;
   p = malloc(1*sizeof(double));
   if (sts!=NULL) {
      *p = value;
   } 

}

static PyMethodDef AccessMethods[] = {
  ...
  {"system", access_system, METH_VARARGS,
   "allocate a double object."},
  ...
  {NULL,NULL,0,NULL} /*Sentinel*/ 

};

static struct PyModuleDef accessmodule = {
  PyModuleDef_HEAD_INIT,
  "access", /*name of the module*/
  access_doc, /* module documentation */
  -1,  /* size of per-interpreter state of the module
         or -1 if the module keeps state in global variables. */
  AccessMethods
};

PyINIT_FUNC
PyInit_access(void) {
   return PyModule_Create(&accessmodule);
}


