#pragma once

#include <cstddef>
#include <iostream>


#define UNIFORM
//#define DEBUGMODE
//#define HARDCORE 


typedef std::size_t SiteIndex;
typedef std::size_t BondIndex;

typedef std::size_t SiteType;
typedef std::pair<SiteType, SiteType> BondType;

const std::size_t ZC=2;
typedef int StateType;

typedef unsigned long index_type;
const double tol = 1e-14;
