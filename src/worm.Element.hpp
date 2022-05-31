/* $Id: worm.hpp,v 1.1 2006/09/09 9:21:44 pollet Exp $ */

#ifndef worm_Element_HPP
#define worm_Element_HPP

#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <list>
#include <stdint.h>

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>



template <typename Site_t, typename State_t, size_t zcoord>
class Element
{
 public:
  Element() {};
  Element(const State_t n0, const State_t n1, const SiteIndex s0, const double t, const int color, const typename std::list< Element<Site_t, State_t, zcoord> >::iterator  v[zcoord] ) {
    mTime = t;
    mBefore = n0;
    mAfter = n1;
    //mFrom = s0;
    //mTo = s1;
    mLink = s0;
    mColor = color;
    for (std::size_t i = 0; i < zcoord; i++) mAssoc[i] = v[i];
  }
  Element(const State_t n0, const State_t n1, const SiteIndex s0, const double t, const int color ) {
    mTime = t;
    mBefore = n0;
    mAfter = n1;
    //mFrom = s0;
    //mTo = s1;
    mLink = s0;
    mColor = color;
    //mAssoc.resize(assoc_size);
  }

  void init(const State_t n0, const State_t n1, const SiteIndex s0, const double t, const int col) {
    mTime = t;
    mBefore = n0;
    mAfter = n1;
    //from(s0);
    //to(s1);
    mLink = s0;
    mColor = col;
  }
  /*
  void init(const State_t n0, const State_t n1, const SiteIndex s0, const double t, const int col, typename std::list<Element<Site_t, State_t, zcoord> >::iterator v[zcoord] ) {
    mTime = t;
    mBefore = n0;
    mAfter = n1;
    //from(s0);
    //to(s1);
    mLink = s0;
    mColor = col;
    for (std::size_t i = 0; i < v.size(); i++) mAssoc[i] = v[i];
  }
  */
  ~Element() {}
  Element(const Element<Site_t, State_t, zcoord>& src) : mTime(src.mTime), mBefore(src.mBefore), mAfter(src.mAfter), mLink(src.mLink),  mColor(src.mColor) {
    for (std::size_t i = 0; i < zcoord; i++) mAssoc[i] = src.mAssoc[i];
  }
  Element& operator=(const Element<Site_t, State_t, zcoord> rhs) {
    if (this == &rhs) return (*this);
    mTime = rhs.mTime;
    mBefore = rhs.mBefore;
    mAfter = rhs.mAfter;
    mLink = rhs.mLink;
    mColor = rhs.mColor;
    for (std::size_t i = 0; i < zcoord; i++) mAssoc[i] = rhs.mAssoc[i];
    return *this;
  }
  friend std::ostream &operator<<( std::ostream &os, const Element& rhs) {
    os  << rhs.color() << "\t"
        << rhs.time() << "\t"
        //<< rhs.from() << "\t"
        //<< rhs.to() << "\t"
        << rhs.link() << "\t"
        << rhs.before() << "\t"
        << rhs.after();
    return (os);
  }
  friend std::istream& operator>>(std::istream& is, Element& e) {
    //is >> e.mColor>> e.mTime >> e.mFrom >> e.mTo >> e.mBefore >> e.mAfter;
    is >> e.mColor>> e.mTime >> e.mLink >> e.mBefore >> e.mAfter;
    return (is);
  }
  friend bool operator== (const Element& lhs, const Element& rhs) {
    return  ( (lhs.mColor == rhs.mColor) && (lhs.mLink == rhs.mLink) && (lhs.mTime == rhs.mTime) && (lhs.mBefore == rhs.mBefore) && (lhs.mAfter == rhs.mAfter) );
  }


  double time() const {return (mTime);}
  State_t before() const {return (mBefore);}
  State_t after() const {return (mAfter);}
  //int name() const { return (mName);}
  int color() const {return (mColor);}
  //SiteType from() const {return mFrom;}
  //SiteType to() const { return mTo;}
  SiteIndex link() const { return mLink;}
  
  void time(const double t) {mTime=t;};
  void before(const State_t n) {mBefore = n;};
  void after(const State_t n) {mAfter = n;};
  //void from(const SiteType s) {mFrom = s;}
  //void to(const SiteType s) { mTo = s;}
  void link(const SiteIndex l)  { mLink = l;}
  //void name(const int n) {mName = n;}
  void color(const int col) {mColor = col;}
  //void assoc(const std::list<Element>::iterator  v[zcoord]) {
  void assoc(const typename std::list<Element<Site_t, State_t, zcoord> >::iterator   v[zcoord]) {
    for (int i = 0; i < zcoord; i++) mAssoc[i] = v[i];
  }

  void print() {
    std::cout << "\n" << this << "\tcolor : " << mColor << "\ttime : " << mTime << "\tlink : " << mLink << "\tbefore : " << mBefore << "\tafter : " << mAfter;
    //for (size_t j=0; j < mAssoc.size(); j++) std::cout << "\t" << &(*mAssoc[j]);
    //for (size_t j=0; j < zcoord; j++) std::cout << "\t" << mAssoc[j]->time() * mAssoc[j]->color() ;
  }
  //void get_assoc(typename std::list<Element<Site_t, State_t, zcoord> >::iterator  v[zcoord]) {
  //  for (std::size_t i = 0; i < mAssoc.size(); i++) {
  //    v[i] = mAssoc[i];
  //  }
  //}
  typename std::list<Element<Site_t, State_t, zcoord> >::iterator& get_assoc(SiteIndex s) { return mAssoc[s];}
  void set_assoc(const SiteIndex j, const typename std::list<Element<Site_t, State_t, zcoord> >::iterator v) {
    mAssoc[j] = v;
  }

  // save to HDF5
  //void save(alps::hdf5::archive& ar) const {
  //    std::vector<T> data(mx, mx+DIM);
  //    ar[""] << data;
  //}

  // load from HDF5
  //void load(alps::hdf5::archive& ar) {
  //    std::vector<T> data;
  //    ar[""] >> data;
  //    std::copy(data.begin(), data.end(), mx);
  //}


 private:
  double mTime;  // time of the interaction
  State_t mBefore;     // occupation on the site before the interaction
  State_t mAfter;      // occupation on the site after the interaction
  Site_t mLink;   // site to which element is linked
  //SiteType mFrom;  // site from which the particle hops
  //SiteType mTo;    // site to which the particle hops
  typename std::list<Element<Site_t, State_t,zcoord> >::iterator mAssoc[zcoord];   // associations; i.e. iterator to the kink on nb sites that are equal or just greater in time
  //int mName;       // numerical name of an element; used to save/load associations of elements to/from file
  int mColor;      // (1) red boson interaction (2) blue boson interaction (-1) red boson worm (-2) blue boson worm (0) dummy etc
};






#endif
