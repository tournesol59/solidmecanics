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
  access_doc

