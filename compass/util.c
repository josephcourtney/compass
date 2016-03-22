#include <Python.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

#define square(x) ((x)*(x))

double mhd_d6( double **a, double **b, int a_len, int b_len, int dim ) {
	double sum = 0.0;
	int i,j,k;
	double dist, min_dist;
	for( i = 0; i < a_len; i++ ) {
        min_dist = DBL_MAX;
		for( j = 0; j < b_len; j++ ) {
            dist = 0.0;
            for( k = 0; k < dim; k++ ) {
                dist += square(a[i][k] - b[j][k]);
            }
			if (dist < min_dist) min_dist = dist;
		}
		sum += sqrt(min_dist);
    }        
    return sum / a_len;
}

double mhd_f2( double a, double b ) {
    return (a > b ? a : b );
}

double mhd( double **a, double **b, int a_len, int b_len, int dim ) {
    return mhd_f2(mhd_d6(a,b,a_len,b_len,dim),mhd_d6(b,a,b_len,a_len,dim));
}

int convert_2d_pylist( PyObject *list, double ***array, int *rows, int *cols ) {
    int i, j;
    PyObject *item;
    PyObject *num;
    // Is the object a list?
    if (!PyList_Check(list))
        return -1;
    *rows = PySequence_Length(list);
    
    // Does the list contain elements?
    if (*rows < 0) {
        return -1;
    }

    item = PySequence_GetItem(list,0);
    if(!PySequence_Check(item)) {
        return -1;
    }
    *cols = PySequence_Length(item);
    if(*cols < 0) {
        return -1;
    }
    Py_DECREF(item);

    *array = (double **)malloc(*rows * sizeof(double*) );
    for( i = 0; i < *rows; i++ ) {
        (*array)[i] = (double*)malloc(*cols*sizeof(double));
    }
        
    
    for( i = 0; i < *rows; i++ ) {
        item = PySequence_GetItem(list,i);
        if(!PySequence_Check(item)) {
            return -1;
        }
        if( PySequence_Length(item) != *cols) {
            return -1;
        }
        for( j = 0; j < *cols; j++ ) {
            num = PySequence_GetItem( item, j );
            if (PyInt_Check(num)) {
                (*array)[i][j] = (double)PyInt_AsLong( num );
            }
            else if (PyLong_Check(num)) {
                (*array)[i][j] = PyLong_AsDouble( num );
            }
            else if (PyFloat_Check(num)) {
                (*array)[i][j] = PyFloat_AsDouble( num );
            }
            else {
                return -1;
            }
        }
        Py_DECREF(item);        
    }
    return 0;
}

static PyObject *calculate_MHD(PyObject *self, PyObject *args) {
    PyObject *list_a, *list_b;
    double **a = NULL;
    double **b = NULL;
    int a_len, b_len, a_dim, b_dim;
    if (!PyArg_ParseTuple(args, "OO", &list_a, &list_b)) {
        printf("Problem Parsing Arguments!\n");
        return NULL;
    }
    if( convert_2d_pylist( list_a, &a, &a_len, &a_dim ) != 0 ) {
        printf("Problem converting list #1\n");
        return NULL;
    }    

    if( convert_2d_pylist( list_b, &b, &b_len, &b_dim) != 0 ) {
        printf("Problem converting list #2\n");
        return NULL;
    }

    if( a_dim != b_dim) {
        printf("Dimension mismatch!\n");
        return NULL;
    }
    return Py_BuildValue("d", mhd( a, b, a_len, b_len, a_dim ) );
}

PyMethodDef compass_util_methods[] = {
    {"mhd", calculate_MHD, METH_VARARGS, "Calculate the Modified Hausdorff Distance between two lists of tuples"},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initutil(void) {
    (void) Py_InitModule("compass.util", compass_util_methods);
}

