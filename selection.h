#ifndef SELECTION_H
#define SELECTION_H

#include "phasespace.h"

#include <iostream>
using std::ostream;

class Selection 
{
public:
  enum Type { GENERATORCUTS, MINIMALCUTS, STANDARDCUTS, HIGGSSEARCH };
  virtual const bool includes(const PhaseSpace& measuredPS, const PhaseSpace& truePS) const = 0;
  virtual const Type type() const = 0;
  virtual void writeSelection(ostream& os) const = 0;
  virtual void writeSelectionCode(ostream& os) const = 0;
  Selection();
  virtual ~Selection() {}
};

// --------------------------------------------------------------------

class MinimalCutsSelection : public Selection
{
public:
  const bool includes(const PhaseSpace& measuredPS, const PhaseSpace& truePS) const;
  const Type type() const;
  void writeSelection(ostream& os) const;
  void writeSelectionCode(ostream& os) const;
  MinimalCutsSelection();
};

inline const Selection::Type MinimalCutsSelection::type() const { return MINIMALCUTS; }

inline void MinimalCutsSelection::writeSelection(ostream& os) const { os << "minimal cuts"; }

// --------------------------------------------------------------------

class GeneratorCutsSelection : public Selection
{
public:
  const bool includes(const PhaseSpace& measuredPS, const PhaseSpace& truePS) const;
    const float NextToLeading(const float pt1, const float pt2, const float pt3, const float pt4) const;
  const Type type() const;
  void writeSelection(ostream& os) const;
  void writeSelectionCode(ostream& os) const;
  GeneratorCutsSelection();
};

inline const Selection::Type GeneratorCutsSelection::type() const { return GENERATORCUTS; }

inline void GeneratorCutsSelection::writeSelection(ostream& os) const { os << "generator cuts"; }

// --------------------------------------------------------------------

class StandardCutsSelection : public Selection
{
public:
  const bool includes(const PhaseSpace& measuredPS, const PhaseSpace& truePS) const;
  const Type type() const;
  void writeSelection(ostream& os) const;
  void writeSelectionCode(ostream& os) const;
  StandardCutsSelection();
};

inline const Selection::Type StandardCutsSelection::type() const { return STANDARDCUTS; }

inline void StandardCutsSelection::writeSelection(ostream& os) const { os << "standard cuts"; }

// --------------------------------------------------------------------

class HiggsSearchSelection : public Selection
{
public:
  const bool includes(const PhaseSpace& measuredPS, const PhaseSpace& truePS) const;
  const Type type() const;
  void writeSelection(ostream& os) const;
  void writeSelectionCode(ostream& os) const;
  HiggsSearchSelection();
};

inline const Selection::Type HiggsSearchSelection::type() const { return HIGGSSEARCH; }

inline void HiggsSearchSelection::writeSelection(ostream& os) const { os << "Higgs search cuts"; }

#endif  /* SELECTION_H */
