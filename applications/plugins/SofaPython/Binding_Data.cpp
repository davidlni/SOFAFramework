/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Plugins                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#include <sofa/core/objectmodel/BaseData.h>
#include <sofa/defaulttype/DataTypeInfo.h>
#include <sofa/core/objectmodel/Data.h>
#include "Binding_LinearSpring.h"

using namespace sofa::core::objectmodel;
using namespace sofa::defaulttype;


// TODO:
// se servir du DataTypeInfo pour utiliser directement les bons type :-)
// Il y a un seul type "Data" exposé en python, le transtypage est géré automatiquement

#include "Binding_Data.h"

extern "C" PyObject * Data_getAttr_name(PyObject *self, void*)
{
    BaseData* data=((PyPtr<BaseData>*)self)->object; // TODO: check dynamic cast
    return PyString_FromString(data->getName().c_str());
}
extern "C" int Data_setAttr_name(PyObject *self, PyObject * args, void*)
{
    BaseData* data=((PyPtr<BaseData>*)self)->object; // TODO: check dynamic cast
    char *str = PyString_AsString(args); // pour les setters, un seul objet et pas un tuple....
    data->setName(str);
    return 0;
}

PyObject *GetDataValuePython(BaseData* data)
{
    // depending on the data type, we return the good python type (int, float, sting, array, ...)

    const AbstractTypeInfo *typeinfo = data->getValueTypeInfo();
    int rowWidth = typeinfo->size();
    int nbRows = typeinfo->size(data->getValueVoidPtr()) / typeinfo->size();

    // special cases...
    Data<sofa::helper::vector<LinearSpring<SReal> > >* vectorLinearSpring = dynamic_cast<Data<sofa::helper::vector<LinearSpring<SReal> > >*>(data);
    if (vectorLinearSpring)
    {
        // special type, a vector of LinearSpring objects
        if (typeinfo->size(data->getValueVoidPtr())==1)
        {
            // this type is NOT a vector; return directly the proper native type
            const LinearSpring<SReal> value = vectorLinearSpring->getValue()[0];
            LinearSpring<SReal> *obj = new LinearSpring<SReal>(value.m1,value.m2,value.ks,value.kd,value.initpos);
            return SP_BUILD_PYPTR(LinearSpring,LinearSpring<SReal>,obj,true); // "true", because I manage the deletion myself
        }
        else
        {
            PyObject *rows = PyList_New(nbRows);
            for (int i=0; i<nbRows; i++)
            {
                PyObject *row = PyList_New(rowWidth);
                for (int j=0; j<rowWidth; j++)
                {
                    // build each value of the list
                    const LinearSpring<SReal> value = vectorLinearSpring->getValue()[i*rowWidth+j];
                    LinearSpring<SReal> *obj = new LinearSpring<SReal>(value.m1,value.m2,value.ks,value.kd,value.initpos);
                    PyList_SetItem(row,j,SP_BUILD_PYPTR(LinearSpring,LinearSpring<SReal>,obj,true));
                }
                PyList_SetItem(rows,i,row);
            }

            return rows;
        }

    }

    if (typeinfo->size(data->getValueVoidPtr())==1)
    {
        // this type is NOT a vector; return directly the proper native type
        if (typeinfo->Text())
        {
            // it's some text
            return PyString_FromString(typeinfo->getTextValue(data->getValueVoidPtr(),0).c_str());
        }
        if (typeinfo->Scalar())
        {
            // it's a SReal
            return PyFloat_FromDouble(typeinfo->getScalarValue(data->getValueVoidPtr(),0));
        }
        if (typeinfo->Integer())
        {
            // it's some Integer...
            return PyInt_FromLong(typeinfo->getIntegerValue(data->getValueVoidPtr(),0));
        }
    }
    else
    {
        // this is a vector; return a python list of the corrsponding type (ints, scalars or strings)

        PyObject *rows = PyList_New(nbRows);
        for (int i=0; i<nbRows; i++)
        {
            PyObject *row = PyList_New(rowWidth);
            for (int j=0; j<rowWidth; j++)
            {
                // build each value of the list
                if (typeinfo->Text())
                {
                    // it's some text
                    PyList_SetItem(row,j,PyString_FromString(typeinfo->getTextValue(data->getValueVoidPtr(),i*rowWidth+j).c_str()));
                }
                else if (typeinfo->Scalar())
                {
                    // it's a SReal
                    PyList_SetItem(row,j,PyFloat_FromDouble(typeinfo->getScalarValue(data->getValueVoidPtr(),i*rowWidth+j)));
                }
                else if (typeinfo->Integer())
                {
                    // it's some Integer...
                    PyList_SetItem(row,j,PyInt_FromLong(typeinfo->getIntegerValue(data->getValueVoidPtr(),i*rowWidth+j)));
                }
                else
                {
                    // this type is not yet supported
                    printf("<SofaPython> BaseData_getAttr_value WARNING: unsupported native type=%s ; returning string value\n",data->getValueTypeString().c_str());
                    PyList_SetItem(row,j,PyString_FromString(typeinfo->getTextValue(data->getValueVoidPtr(),i*rowWidth+j).c_str()));
                }
            }
            PyList_SetItem(rows,i,row);
        }

        return rows;
    }
    // default (should not happen)...
    printf("<SofaPython> BaseData_getAttr_value WARNING: unsupported native type=%s ; returning string value\n",data->getValueTypeString().c_str());
    return PyString_FromString(data->getValueString().c_str());
}

bool SetDataValuePython(BaseData* data, PyObject* args)
{
    // de quel type est args ?
    bool isInt = PyInt_Check(args);
    bool isScalar = PyFloat_Check(args);
    bool isString = PyString_Check(args);
    bool isList = PyList_Check(args);
    const AbstractTypeInfo *typeinfo = data->getValueTypeInfo(); // info about the data value
    int rowWidth = typeinfo->size();
    int nbRows = typeinfo->size(data->getValueVoidPtr()) / typeinfo->size();

    // special cases...
    Data<sofa::helper::vector<LinearSpring<SReal> > >* dataVectorLinearSpring = dynamic_cast<Data<sofa::helper::vector<LinearSpring<SReal> > >*>(data);
    if (dataVectorLinearSpring)
    {
        // special type, a vector of LinearSpring objects

        if (!isList)
        {
            // one value
            // check the python object type
            if (rowWidth*nbRows<1 || !PyObject_IsInstance(args,reinterpret_cast<PyObject*>(&SP_SOFAPYTYPEOBJECT(LinearSpring))))
            {
                // type mismatch or too long list
                PyErr_BadArgument();
                return false;
            }

            LinearSpring<SReal>* obj=dynamic_cast<LinearSpring<SReal>*>(((PyPtr<LinearSpring<SReal> >*)args)->object);
            sofa::helper::vector<LinearSpring<SReal> >* vectorLinearSpring = dataVectorLinearSpring->beginEdit();

            (*vectorLinearSpring)[0].m1 = obj->m1;
            (*vectorLinearSpring)[0].m2 = obj->m2;
            (*vectorLinearSpring)[0].ks = obj->ks;
            (*vectorLinearSpring)[0].kd = obj->kd;
            (*vectorLinearSpring)[0].initpos = obj->initpos;

            dataVectorLinearSpring->endEdit();

            return true;
        }
        else
        {
            // values list
            // is-it a double-dimension list ?
            //PyObject *firstRow = PyList_GetItem(args,0);

            if (PyList_Check(PyList_GetItem(args,0)))
            {
                // two-dimension array!

                // right number if rows ?
                if (PyList_Size(args)!=nbRows)
                {
                    // only a warning; do not raise an exception...
                    printf("<SofaPython> Warning: list size mismatch for data \"%s\" (incorrect rows count)\n",data->getName().c_str());
                    if (PyList_Size(args)<nbRows)
                        nbRows = PyList_Size(args);
                }

                sofa::helper::vector<LinearSpring<SReal> >* vectorLinearSpring = dataVectorLinearSpring->beginEdit();

                // let's fill our rows!
                for (int i=0; i<nbRows; i++)
                {
                    PyObject *row = PyList_GetItem(args,i);

                    // right number if list members ?
                    int size = rowWidth;
                    if (PyList_Size(row)!=size)
                    {
                        // only a warning; do not raise an exception...
                        printf("<SofaPython> Warning: row %i size mismatch for data \"%s\"\n",i,data->getName().c_str());
                        if (PyList_Size(row)<size)
                            size = PyList_Size(row);
                    }

                    // okay, let's set our list...
                    for (int j=0; j<size; j++)
                    {

                        PyObject *listElt = PyList_GetItem(row,j);
                        if(!PyObject_IsInstance(listElt,reinterpret_cast<PyObject*>(&SP_SOFAPYTYPEOBJECT(LinearSpring))))
                        {
                            // type mismatch
                            dataVectorLinearSpring->endEdit();
                            PyErr_BadArgument();
                            return false;
                        }
                        LinearSpring<SReal>* spring=dynamic_cast<LinearSpring<SReal>*>(((PyPtr<LinearSpring<SReal> >*)listElt)->object);


                        (*vectorLinearSpring)[j+i*rowWidth].m1 = spring->m1;
                        (*vectorLinearSpring)[j+i*rowWidth].m2 = spring->m2;
                        (*vectorLinearSpring)[j+i*rowWidth].ks = spring->ks;
                        (*vectorLinearSpring)[j+i*rowWidth].kd = spring->kd;
                        (*vectorLinearSpring)[j+i*rowWidth].initpos = spring->initpos;

                    }



                }

                dataVectorLinearSpring->endEdit();

                return true;

            }
            else
            {
                // it is a one-dimension only array
                // right number if list members ?
                int size = rowWidth*nbRows;
                if (PyList_Size(args)!=size)
                {
                    // only a warning; do not raise an exception...
                    printf("<SofaPython> Warning: list size mismatch for data \"%s\"\n",data->getName().c_str());
                    if (PyList_Size(args)<size)
                        size = PyList_Size(args);
                }

                sofa::helper::vector<LinearSpring<SReal> >* vectorLinearSpring = dataVectorLinearSpring->beginEdit();

                // okay, let's set our list...
                for (int i=0; i<size; i++)
                {

                    PyObject *listElt = PyList_GetItem(args,i);

                    if(!PyObject_IsInstance(listElt,reinterpret_cast<PyObject*>(&SP_SOFAPYTYPEOBJECT(LinearSpring))))
                    {
                        // type mismatch
                        dataVectorLinearSpring->endEdit();
                        PyErr_BadArgument();
                        return false;
                    }

                    LinearSpring<SReal>* spring=dynamic_cast<LinearSpring<SReal>*>(((PyPtr<LinearSpring<SReal> >*)listElt)->object);

                    (*vectorLinearSpring)[i].m1 = spring->m1;
                    (*vectorLinearSpring)[i].m2 = spring->m2;
                    (*vectorLinearSpring)[i].ks = spring->ks;
                    (*vectorLinearSpring)[i].kd = spring->kd;
                    (*vectorLinearSpring)[i].initpos = spring->initpos;


    /*
                    if (PyFloat_Check(listElt))
                    {
                        // it's a scalar
                        if (!typeinfo->Scalar())
                        {
                            // type mismatch
                            PyErr_BadArgument();
                            return false;
                        }
                        SReal value = PyFloat_AsDouble(listElt);
                        typeinfo->setScalarValue((void*)data->getValueVoidPtr(),i,value);
                    }
     */
                }
                dataVectorLinearSpring->endEdit();

                return true;
            }
        }


        return false;
    }


    if (isInt)
    {
        // it's an int

        if (rowWidth*nbRows<1 || (!typeinfo->Integer() && !typeinfo->Scalar()))
        {
            // type mismatch or too long list
            PyErr_BadArgument();
            return false;
        }
        long value = PyInt_AsLong(args);
        if (typeinfo->Scalar())
            typeinfo->setScalarValue((void*)data->getValueVoidPtr(),0,(SReal)value); // cast int to float
        else
            typeinfo->setIntegerValue((void*)data->getValueVoidPtr(),0,value);
        return true;
    }
    else if (isScalar)
    {
        // it's a scalar
        if (rowWidth*nbRows<1 || !typeinfo->Scalar())
        {
            // type mismatch or too long list
            PyErr_BadArgument();
            return false;
        }
        SReal value = PyFloat_AsDouble(args);
        typeinfo->setScalarValue((void*)data->getValueVoidPtr(),0,value);
        return true;
    }
    else if (isString)
    {
        // it's a string
        char *str = PyString_AsString(args); // pour les setters, un seul objet et pas un tuple....
        data->read(str);
        return true;
    }
    else if (isList)
    {
        // it's a list
        // check list emptyness
        if (PyList_Size(args)==0)
        {
            // empty list: ignored
            return true;
        }

        // is-it a double-dimension list ?
        //PyObject *firstRow = PyList_GetItem(args,0);

        if (PyList_Check(PyList_GetItem(args,0)))
        {
            // two-dimension array!

            // right number if rows ?
            if (PyList_Size(args)!=nbRows)
            {
                // only a warning; do not raise an exception...
                printf("<SofaPython> Warning: list size mismatch for data \"%s\" (incorrect rows count)\n",data->getName().c_str());
                if (PyList_Size(args)<nbRows)
                    nbRows = PyList_Size(args);
            }

            // let's fill our rows!
            for (int i=0; i<nbRows; i++)
            {
                PyObject *row = PyList_GetItem(args,i);

                // right number if list members ?
                int size = rowWidth;
                if (PyList_Size(row)!=size)
                {
                    // only a warning; do not raise an exception...
                    printf("<SofaPython> Warning: row %i size mismatch for data \"%s\"\n",i,data->getName().c_str());
                    if (PyList_Size(row)<size)
                        size = PyList_Size(row);
                }

                // okay, let's set our list...
                for (int j=0; j<size; j++)
                {

                    PyObject *listElt = PyList_GetItem(row,j);

                    if (PyInt_Check(listElt))
                    {
                        // it's an int
                        if (typeinfo->Integer())
                        {
                            // integer value
                            long value = PyInt_AsLong(listElt);
                            typeinfo->setIntegerValue((void*)data->getValueVoidPtr(),i*rowWidth+j,value);
                        }
                        else if (typeinfo->Scalar())
                        {
                            // cast to scalar value
                            SReal value = (SReal)PyInt_AsLong(listElt);
                            typeinfo->setScalarValue((void*)data->getValueVoidPtr(),i*rowWidth+j,value);
                        }
                        else
                        {
                            // type mismatch
                            PyErr_BadArgument();
                            return false;
                        }
                    }
                    else if (PyFloat_Check(listElt))
                    {
                        // it's a scalar
                        if (!typeinfo->Scalar())
                        {
                            // type mismatch
                            PyErr_BadArgument();
                            return false;
                        }
                        SReal value = PyFloat_AsDouble(listElt);
                        typeinfo->setScalarValue((void*)data->getValueVoidPtr(),i*rowWidth+j,value);
                    }
                    else if (PyString_Check(listElt))
                    {
                        // it's a string
                        if (!typeinfo->Text())
                        {
                            // type mismatch
                            PyErr_BadArgument();
                            return false;
                        }
                        char *str = PyString_AsString(listElt); // pour les setters, un seul objet et pas un tuple....
                        typeinfo->setTextValue((void*)data->getValueVoidPtr(),i*rowWidth+j,str);
                    }
                    else
                    {
                        printf("Lists not yet supported...\n");
                        PyErr_BadArgument();
                        return false;

                    }
                }



            }
            return true;

        }
        else
        {
            // it is a one-dimension only array
            // right number if list members ?
            int size = rowWidth*nbRows;
            if (PyList_Size(args)!=size)
            {
                // only a warning; do not raise an exception...
                printf("<SofaPython> Warning: list size mismatch for data \"%s\"\n",data->getName().c_str());
                if (PyList_Size(args)<size)
                    size = PyList_Size(args);
            }

            // okay, let's set our list...
            for (int i=0; i<size; i++)
            {

                PyObject *listElt = PyList_GetItem(args,i);

                if (PyInt_Check(listElt))
                {
                    // it's an int
                    if (typeinfo->Integer())
                    {
                        // integer value
                        long value = PyInt_AsLong(listElt);
                        typeinfo->setIntegerValue((void*)data->getValueVoidPtr(),i,value);
                    }
                    else if (typeinfo->Scalar())
                    {
                        // cast to scalar value
                        SReal value = (SReal)PyInt_AsLong(listElt);
                        typeinfo->setScalarValue((void*)data->getValueVoidPtr(),i,value);
                    }
                    else
                    {
                        // type mismatch
                        PyErr_BadArgument();
                        return false;
                    }
                }
                else if (PyFloat_Check(listElt))
                {
                    // it's a scalar
                    if (!typeinfo->Scalar())
                    {
                        // type mismatch
                        PyErr_BadArgument();
                        return false;
                    }
                    SReal value = PyFloat_AsDouble(listElt);
                    typeinfo->setScalarValue((void*)data->getValueVoidPtr(),i,value);
                }
                else if (PyString_Check(listElt))
                {
                    // it's a string
                    if (!typeinfo->Text())
                    {
                        // type mismatch
                        PyErr_BadArgument();
                        return false;
                    }
                    char *str = PyString_AsString(listElt); // pour les setters, un seul objet et pas un tuple....
                    typeinfo->setTextValue((void*)data->getValueVoidPtr(),i,str);
                }
                else
                {
                    printf("Lists not yet supported...\n");
                    PyErr_BadArgument();
                    return false;

                }
            }

            return true;
        }

    }

    return false;

}



extern "C" PyObject * Data_getAttr_value(PyObject *self, void*)
{
    BaseData* data=((PyPtr<BaseData>*)self)->object; // TODO: check dynamic cast
    return GetDataValuePython(data);
}

extern "C" int Data_setAttr_value(PyObject *self, PyObject * args, void*)
{
    BaseData* data=((PyPtr<BaseData>*)self)->object; // TODO: check dynamic cast
    if (SetDataValuePython(data,args))
        return 0;   // OK


    printf("<SofaPython> argument type not supported\n");
    PyErr_BadArgument();
    return -1;
}

// access ONE element of the vector
extern "C" PyObject * Data_getValue(PyObject *self, PyObject * args)
{
    BaseData* data=((PyPtr<BaseData>*)self)->object;
    const AbstractTypeInfo *typeinfo = data->getValueTypeInfo(); // info about the data value
    int index;
    if (!PyArg_ParseTuple(args, "i",&index))
    {
        PyErr_BadArgument();
        return 0;
    }
    if ((unsigned int)index>=typeinfo->size())
    {
        // out of bounds!
        printf("<SofaPython> Error: Data.getValue index overflow\n");
        PyErr_BadArgument();
        return 0;
    }
    if (typeinfo->Scalar())
        return PyFloat_FromDouble(typeinfo->getScalarValue(data->getValueVoidPtr(),index));
    if (typeinfo->Integer())
        return PyInt_FromLong(typeinfo->getIntegerValue(data->getValueVoidPtr(),index));
    if (typeinfo->Text())
        return PyString_FromString(typeinfo->getTextValue(data->getValueVoidPtr(),index).c_str());

    // should never happen....
    printf("<SofaPython> Error: Data.getValue unknown data type\n");
    PyErr_BadArgument();
    return 0;
}
extern "C" PyObject * Data_setValue(PyObject *self, PyObject * args)
{
    BaseData* data=((PyPtr<BaseData>*)self)->object;
    const AbstractTypeInfo *typeinfo = data->getValueTypeInfo(); // info about the data value
    int index;
    PyObject *value;
    if (!PyArg_ParseTuple(args, "iO",&index,&value))
    {
        PyErr_BadArgument();
        return 0;
    }
    if ((unsigned int)index>=typeinfo->size())
    {
        // out of bounds!
        printf("<SofaPython> Error: Data.setValue index overflow\n");
        PyErr_BadArgument();
        return 0;
    }
    if (typeinfo->Scalar() && PyFloat_Check(value))
    {
        typeinfo->setScalarValue((void*)data->getValueVoidPtr(),index,PyFloat_AsDouble(value));
        return PyInt_FromLong(0);
    }
    if (typeinfo->Integer() && PyInt_Check(value))
    {
        typeinfo->setIntegerValue((void*)data->getValueVoidPtr(),index,PyInt_AsLong(value));
        return PyInt_FromLong(0);
    }
    if (typeinfo->Text() && PyString_Check(value))
    {
        typeinfo->setTextValue((void*)data->getValueVoidPtr(),index,PyString_AsString(value));
        return PyInt_FromLong(0);
    }

    // should never happen....
    printf("<SofaPython> Error: Data.setValue type mismatch\n");
    PyErr_BadArgument();
    return 0;
}


extern "C" PyObject * Data_getValueTypeString(PyObject *self, PyObject * /*args*/)
{
    BaseData* data=((PyPtr<BaseData>*)self)->object;
    return PyString_FromString(data->getValueTypeString().c_str());
}

extern "C" PyObject * Data_getValueString(PyObject *self, PyObject * /*args*/)
{
    BaseData* data=((PyPtr<BaseData>*)self)->object;
    return PyString_FromString(data->getValueString().c_str());
}

SP_CLASS_METHODS_BEGIN(Data)
SP_CLASS_METHOD(Data,getValueTypeString)
SP_CLASS_METHOD(Data,getValueString)
SP_CLASS_METHOD(Data,setValue)
SP_CLASS_METHOD(Data,getValue)
SP_CLASS_METHODS_END


SP_CLASS_ATTRS_BEGIN(Data)
SP_CLASS_ATTR(Data,name)
//SP_CLASS_ATTR(BaseData,owner)
SP_CLASS_ATTR(Data,value)
SP_CLASS_ATTRS_END

SP_CLASS_TYPE_BASE_PTR_ATTR(Data,BaseData)
