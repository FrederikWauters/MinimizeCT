// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME DataClassesDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
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

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "TData.h"
#include "DataManager.h"
#include "Globals.h"
#include "Analyzer.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TData(void *p = 0);
   static void *newArray_TData(Long_t size, void *p);
   static void delete_TData(void *p);
   static void deleteArray_TData(void *p);
   static void destruct_TData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TData*)
   {
      ::TData *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TData >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TData", ::TData::Class_Version(), "TData.h", 12,
                  typeid(::TData), DefineBehavior(ptr, ptr),
                  &::TData::Dictionary, isa_proxy, 4,
                  sizeof(::TData) );
      instance.SetNew(&new_TData);
      instance.SetNewArray(&newArray_TData);
      instance.SetDelete(&delete_TData);
      instance.SetDeleteArray(&deleteArray_TData);
      instance.SetDestructor(&destruct_TData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TData*)
   {
      return GenerateInitInstanceLocal((::TData*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TData*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TData::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TData::Class_Name()
{
   return "TData";
}

//______________________________________________________________________________
const char *TData::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TData*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TData::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TData*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TData::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TData*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TData::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TData*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TData::Streamer(TBuffer &R__b)
{
   // Stream an object of class TData.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TData::Class(),this);
   } else {
      R__b.WriteClassBuffer(TData::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TData(void *p) {
      return  p ? new(p) ::TData : new ::TData;
   }
   static void *newArray_TData(Long_t nElements, void *p) {
      return p ? new(p) ::TData[nElements] : new ::TData[nElements];
   }
   // Wrapper around operator delete
   static void delete_TData(void *p) {
      delete ((::TData*)p);
   }
   static void deleteArray_TData(void *p) {
      delete [] ((::TData*)p);
   }
   static void destruct_TData(void *p) {
      typedef ::TData current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TData

namespace {
  void TriggerDictionaryInitialization_DataClassesDict_Impl() {
    static const char* headers[] = {
"TData.h",
"DataManager.h",
"Globals.h",
"Analyzer.h",
0
    };
    static const char* includePaths[] = {
"/home/frederik/root/include",
"/home/frederik/Analysis/CTMin/New/MinimizeCT/",
0
    };
    static const char* fwdDeclCode = 
"class TData;"
;
    static const char* payloadCode = 
"\n"
"#ifndef G__VECTOR_HAS_CLASS_ITERATOR\n"
"  #define G__VECTOR_HAS_CLASS_ITERATOR\n"
"#endif\n"
"\n"
"#define _BACKWARD_BACKWARD_WARNING_H\n"
"#include \"TData.h\"\n"
"#include \"DataManager.h\"\n"
"#include \"Globals.h\"\n"
"#include \"Analyzer.h\"\n"
"\n"
"#undef  _BACKWARD_BACKWARD_WARNING_H\n"
;
    static const char* classesHeaders[]={
"TData", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("DataClassesDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_DataClassesDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_DataClassesDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_DataClassesDict() {
  TriggerDictionaryInitialization_DataClassesDict_Impl();
}
