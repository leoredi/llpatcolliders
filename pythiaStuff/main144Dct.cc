// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME main144Dct
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "main144Dct.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_RootParticle(void *p = nullptr);
   static void *newArray_RootParticle(Long_t size, void *p);
   static void delete_RootParticle(void *p);
   static void deleteArray_RootParticle(void *p);
   static void destruct_RootParticle(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RootParticle*)
   {
      ::RootParticle *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RootParticle >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("RootParticle", ::RootParticle::Class_Version(), "main144Dct.h", 18,
                  typeid(::RootParticle), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RootParticle::Dictionary, isa_proxy, 4,
                  sizeof(::RootParticle) );
      instance.SetNew(&new_RootParticle);
      instance.SetNewArray(&newArray_RootParticle);
      instance.SetDelete(&delete_RootParticle);
      instance.SetDeleteArray(&deleteArray_RootParticle);
      instance.SetDestructor(&destruct_RootParticle);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RootParticle*)
   {
      return GenerateInitInstanceLocal(static_cast<::RootParticle*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::RootParticle*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_RootEvent(void *p = nullptr);
   static void *newArray_RootEvent(Long_t size, void *p);
   static void delete_RootEvent(void *p);
   static void deleteArray_RootEvent(void *p);
   static void destruct_RootEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RootEvent*)
   {
      ::RootEvent *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RootEvent >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("RootEvent", ::RootEvent::Class_Version(), "main144Dct.h", 45,
                  typeid(::RootEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RootEvent::Dictionary, isa_proxy, 4,
                  sizeof(::RootEvent) );
      instance.SetNew(&new_RootEvent);
      instance.SetNewArray(&newArray_RootEvent);
      instance.SetDelete(&delete_RootEvent);
      instance.SetDeleteArray(&deleteArray_RootEvent);
      instance.SetDestructor(&destruct_RootEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RootEvent*)
   {
      return GenerateInitInstanceLocal(static_cast<::RootEvent*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::RootEvent*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr RootParticle::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *RootParticle::Class_Name()
{
   return "RootParticle";
}

//______________________________________________________________________________
const char *RootParticle::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RootParticle*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int RootParticle::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RootParticle*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RootParticle::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RootParticle*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RootParticle::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RootParticle*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr RootEvent::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *RootEvent::Class_Name()
{
   return "RootEvent";
}

//______________________________________________________________________________
const char *RootEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RootEvent*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int RootEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RootEvent*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RootEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RootEvent*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RootEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RootEvent*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void RootParticle::Streamer(TBuffer &R__b)
{
   // Stream an object of class RootParticle.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RootParticle::Class(),this);
   } else {
      R__b.WriteClassBuffer(RootParticle::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_RootParticle(void *p) {
      return  p ? new(p) ::RootParticle : new ::RootParticle;
   }
   static void *newArray_RootParticle(Long_t nElements, void *p) {
      return p ? new(p) ::RootParticle[nElements] : new ::RootParticle[nElements];
   }
   // Wrapper around operator delete
   static void delete_RootParticle(void *p) {
      delete (static_cast<::RootParticle*>(p));
   }
   static void deleteArray_RootParticle(void *p) {
      delete [] (static_cast<::RootParticle*>(p));
   }
   static void destruct_RootParticle(void *p) {
      typedef ::RootParticle current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::RootParticle

//______________________________________________________________________________
void RootEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class RootEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RootEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(RootEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_RootEvent(void *p) {
      return  p ? new(p) ::RootEvent : new ::RootEvent;
   }
   static void *newArray_RootEvent(Long_t nElements, void *p) {
      return p ? new(p) ::RootEvent[nElements] : new ::RootEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_RootEvent(void *p) {
      delete (static_cast<::RootEvent*>(p));
   }
   static void deleteArray_RootEvent(void *p) {
      delete [] (static_cast<::RootEvent*>(p));
   }
   static void destruct_RootEvent(void *p) {
      typedef ::RootEvent current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::RootEvent

namespace ROOT {
   static TClass *vectorlERootParticlegR_Dictionary();
   static void vectorlERootParticlegR_TClassManip(TClass*);
   static void *new_vectorlERootParticlegR(void *p = nullptr);
   static void *newArray_vectorlERootParticlegR(Long_t size, void *p);
   static void delete_vectorlERootParticlegR(void *p);
   static void deleteArray_vectorlERootParticlegR(void *p);
   static void destruct_vectorlERootParticlegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RootParticle>*)
   {
      vector<RootParticle> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RootParticle>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RootParticle>", -2, "vector", 389,
                  typeid(vector<RootParticle>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERootParticlegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<RootParticle>) );
      instance.SetNew(&new_vectorlERootParticlegR);
      instance.SetNewArray(&newArray_vectorlERootParticlegR);
      instance.SetDelete(&delete_vectorlERootParticlegR);
      instance.SetDeleteArray(&deleteArray_vectorlERootParticlegR);
      instance.SetDestructor(&destruct_vectorlERootParticlegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RootParticle> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<RootParticle>","std::__1::vector<RootParticle, std::__1::allocator<RootParticle>>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<RootParticle>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERootParticlegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<RootParticle>*>(nullptr))->GetClass();
      vectorlERootParticlegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERootParticlegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERootParticlegR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<RootParticle> : new vector<RootParticle>;
   }
   static void *newArray_vectorlERootParticlegR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<RootParticle>[nElements] : new vector<RootParticle>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERootParticlegR(void *p) {
      delete (static_cast<vector<RootParticle>*>(p));
   }
   static void deleteArray_vectorlERootParticlegR(void *p) {
      delete [] (static_cast<vector<RootParticle>*>(p));
   }
   static void destruct_vectorlERootParticlegR(void *p) {
      typedef vector<RootParticle> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<RootParticle>

namespace ROOT {
   // Registration Schema evolution read functions
   int RecordReadRules_main144Dct() {
      return 0;
   }
   static int _R__UNIQUE_DICT_(ReadRules_main144Dct) = RecordReadRules_main144Dct();R__UseDummy(_R__UNIQUE_DICT_(ReadRules_main144Dct));
} // namespace ROOT
namespace {
  void TriggerDictionaryInitialization_main144Dct_Impl() {
    static const char* headers[] = {
"main144Dct.h",
nullptr
    };
    static const char* includePaths[] = {
"../../pythia8315/include",
"/opt/homebrew/Cellar/root/6.36.04/include/root",
"/Users/fredi/cernbox/Physics/llpatcolliders/llpatcolliders/pythiaStuff/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "main144Dct dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$main144Dct.h")))  RootParticle;
class __attribute__((annotate("$clingAutoload$main144Dct.h")))  RootEvent;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "main144Dct dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "main144Dct.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"RootEvent", payloadCode, "@",
"RootParticle", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("main144Dct",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_main144Dct_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_main144Dct_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_main144Dct() {
  TriggerDictionaryInitialization_main144Dct_Impl();
}
