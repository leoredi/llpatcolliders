#ifndef TEST_H
#define TEST_H

#include "TObject.h"

class MyClass : public TObject {
public:
    MyClass() {}
    virtual ~MyClass() {}

    ClassDef(MyClass, 1);
};

#endif // TEST_H
