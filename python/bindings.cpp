#include "vcfx_core.h"
#include <Python.h>
#include <sstream>
#include <string>
#include <vector>

// Helper to convert std::vector<std::string> to Python list
static PyObject *to_py_list(const std::vector<std::string> &vec) {
    PyObject *list = PyList_New(vec.size());
    if (!list)
        return nullptr;
    for (size_t i = 0; i < vec.size(); ++i) {
        PyObject *item = PyUnicode_FromString(vec[i].c_str());
        if (!item) {
            Py_DECREF(list);
            return nullptr;
        }
        PyList_SET_ITEM(list, i, item); // steals reference
    }
    return list;
}

static PyObject *py_trim(PyObject *, PyObject *args) {
    const char *text;
    if (!PyArg_ParseTuple(args, "s", &text))
        return nullptr;
    std::string result = vcfx::trim(text);
    return PyUnicode_FromString(result.c_str());
}

static PyObject *py_split(PyObject *, PyObject *args) {
    const char *text;
    const char *delim;
    if (!PyArg_ParseTuple(args, "ss", &text, &delim))
        return nullptr;
    std::vector<std::string> parts = vcfx::split(text, delim[0]);
    return to_py_list(parts);
}

static PyObject *py_read_file(PyObject *, PyObject *args) {
    const char *path;
    if (!PyArg_ParseTuple(args, "s", &path))
        return nullptr;
    std::string out;
    if (!vcfx::read_file_maybe_compressed(path, out)) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to read file");
        return nullptr;
    }
    return PyBytes_FromStringAndSize(out.data(), out.size());
}

static PyObject *py_get_version(PyObject *, PyObject *) {
    std::string ver = vcfx::get_version();
    return PyUnicode_FromString(ver.c_str());
}

static PyObject *py_read_stream(PyObject *, PyObject *args) {
    Py_buffer buf;
    if (!PyArg_ParseTuple(args, "y*", &buf))
        return nullptr;
    std::string data(static_cast<const char *>(buf.buf), buf.len);
    PyBuffer_Release(&buf);
    std::istringstream ss(data);
    std::string out;
    if (!vcfx::read_maybe_compressed(ss, out)) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to read data");
        return nullptr;
    }
    return PyBytes_FromStringAndSize(out.data(), out.size());
}

static PyMethodDef VcfxMethods[] = {
    {"trim", py_trim, METH_VARARGS, "Trim leading and trailing whitespace"},
    {"split", py_split, METH_VARARGS, "Split a string on the given delimiter"},
    {"read_file_maybe_compressed", py_read_file, METH_VARARGS,
     "Read a (possibly compressed) file and return its contents"},
    {"read_maybe_compressed", py_read_stream, METH_VARARGS, "Decompress bytes if needed and return the contents"},
    {"get_version", py_get_version, METH_NOARGS, "Return VCFX version string"},
    {nullptr, nullptr, 0, nullptr}};

static struct PyModuleDef moduledef = {PyModuleDef_HEAD_INIT, "_vcfx", "Python bindings for VCFX helper functions", -1,
                                       VcfxMethods};

PyMODINIT_FUNC PyInit__vcfx(void) { return PyModule_Create(&moduledef); }
