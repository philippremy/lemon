/* -*- mode: C++; indent-tabs-mode: nil; -*-
 *
 * This file is a part of LEMON, a generic C++ optimization library.
 *
 * Copyright (C) 2003-2013
 * Egervary Jeno Kombinatorikus Optimalizalasi Kutatocsoport
 * (Egervary Research Group on Combinatorial Optimization, EGRES).
 *
 * Permission to use, modify and distribute this software is granted
 * provided that this copyright notice appears in all copies. For
 * precise terms see the accompanying LICENSE file.
 *
 * This software is provided "AS IS" with no warranty of any kind,
 * express or implied, and with no claim as to its suitability for any
 * purpose.
 *
 */

#include <iostream>
#include <lemon/soplex.h>

///\file
///\brief Implementation of the LEMON-SOPLEX lp solver interface.
namespace lemon {

  SoplexLp::SoplexLp() {
    soplex = new soplex::SoPlex;
    messageLevel(MESSAGE_NOTHING);
  }

  SoplexLp::~SoplexLp() {
    delete soplex;
  }

  SoplexLp::SoplexLp(const SoplexLp& lp) {
    rows = lp.rows;
    cols = lp.cols;

    soplex = new soplex::SoPlex;

    _col_names = lp._col_names;
    _col_names_ref = lp._col_names_ref;

    _row_names = lp._row_names;
    _row_names_ref = lp._row_names_ref;

    messageLevel(MESSAGE_NOTHING);
  }

  void SoplexLp::_clear_temporals() {
    _primal_values.clear();
    _dual_values.clear();
  }

  SoplexLp* SoplexLp::newSolver() const {
    SoplexLp* newlp = new SoplexLp();
    return newlp;
  }

  SoplexLp* SoplexLp::cloneSolver() const {
    SoplexLp* newlp = new SoplexLp(*this);
    return newlp;
  }

  const char* SoplexLp::_solverName() const { return "SoplexLp"; }

  int SoplexLp::_addCol() {
    soplex::LPCol c;
    c.setLower(-soplex::infinity);
    c.setUpper(soplex::infinity);
    soplex->addColReal(c);

    _col_names.push_back(std::string());

    return soplex->numColsReal() - 1;
  }

  int SoplexLp::_addRow() {
    soplex::LPRow r;
    r.setLhs(-soplex::infinity);
    r.setRhs(soplex::infinity);
    soplex->addRowReal(r);

    _row_names.push_back(std::string());

    return soplex->numRowsReal() - 1;
  }

  int SoplexLp::_addRow(Value l, ExprIterator b, ExprIterator e, Value u) {
    soplex::DSVector v;
    for (ExprIterator it = b; it != e; ++it) {
      v.add(it->first, it->second);
    }
    soplex::LPRow r(l, v, u);
    soplex->addRowReal(r);

    _row_names.push_back(std::string());

    return soplex->numRowsReal() - 1;
  }


  void SoplexLp::_eraseCol(int i) {
    soplex->removeColReal(i);
    _col_names_ref.erase(_col_names[i]);
    _col_names[i] = _col_names.back();
    _col_names_ref[_col_names.back()] = i;
    _col_names.pop_back();
  }

  void SoplexLp::_eraseRow(int i) {
    soplex->removeRowReal(i);
    _row_names_ref.erase(_row_names[i]);
    _row_names[i] = _row_names.back();
    _row_names_ref[_row_names.back()] = i;
    _row_names.pop_back();
  }

  void SoplexLp::_eraseColId(int i) {
    cols.eraseIndex(i);
    cols.relocateIndex(i, cols.maxIndex());
  }
  void SoplexLp::_eraseRowId(int i) {
    rows.eraseIndex(i);
    rows.relocateIndex(i, rows.maxIndex());
  }

  void SoplexLp::_getColName(int c, std::string &name) const {
    name = _col_names[c];
  }

  void SoplexLp::_setColName(int c, const std::string &name) {
    _col_names_ref.erase(_col_names[c]);
    _col_names[c] = name;
    if (!name.empty()) {
      _col_names_ref.insert(std::make_pair(name, c));
    }
  }

  int SoplexLp::_colByName(const std::string& name) const {
    std::map<std::string, int>::const_iterator it =
      _col_names_ref.find(name);
    if (it != _col_names_ref.end()) {
      return it->second;
    } else {
      return -1;
    }
  }

  void SoplexLp::_getRowName(int r, std::string &name) const {
    name = _row_names[r];
  }

  void SoplexLp::_setRowName(int r, const std::string &name) {
    _row_names_ref.erase(_row_names[r]);
    _row_names[r] = name;
    if (!name.empty()) {
      _row_names_ref.insert(std::make_pair(name, r));
    }
  }

  int SoplexLp::_rowByName(const std::string& name) const {
    std::map<std::string, int>::const_iterator it =
      _row_names_ref.find(name);
    if (it != _row_names_ref.end()) {
      return it->second;
    } else {
      return -1;
    }
  }


  void SoplexLp::_setRowCoeffs(int i, ExprIterator b, ExprIterator e) {
    for (int j = 0; j < soplex->numColsReal(); ++j) {
      soplex->changeElementReal(i, j, 0.0);
    }
    for(ExprIterator it = b; it != e; ++it) {
      soplex->changeElementReal(i, it->first, it->second);
    }
  }

  void SoplexLp::_getRowCoeffs(int i, InsertIterator b) const {
    const soplex::SVector& vec = soplex->rowVectorRealInternal(i);
    for (int k = 0; k < vec.size(); ++k) {
      *b = std::make_pair(vec.index(k), vec.value(k));
      ++b;
    }
  }

  void SoplexLp::_setColCoeffs(int j, ExprIterator b, ExprIterator e) {
    for (int i = 0; i < soplex->numRowsReal(); ++i) {
      soplex->changeElementReal(i, j, 0.0);
    }
    for(ExprIterator it = b; it != e; ++it) {
      soplex->changeElementReal(it->first, j, it->second);
    }
  }

  void SoplexLp::_getColCoeffs(int i, InsertIterator b) const {
    const soplex::SVector& vec = soplex->colVectorRealInternal(i);
    for (int k = 0; k < vec.size(); ++k) {
      *b = std::make_pair(vec.index(k), vec.value(k));
      ++b;
    }
  }

  void SoplexLp::_setCoeff(int i, int j, Value value) {
    soplex->changeElementReal(i, j, value);
  }

  SoplexLp::Value SoplexLp::_getCoeff(int i, int j) const {
    return soplex->rowVectorRealInternal(i)[j];
  }

  void SoplexLp::_setColLowerBound(int i, Value value) {
    LEMON_ASSERT(value != INF, "Invalid bound");
    soplex->changeLowerReal(i, value != -INF ? value : -soplex::infinity);
  }

  SoplexLp::Value SoplexLp::_getColLowerBound(int i) const {
    double value = soplex->lowerReal(i);
    return value != -soplex::infinity ? value : -INF;
  }

  void SoplexLp::_setColUpperBound(int i, Value value) {
    LEMON_ASSERT(value != -INF, "Invalid bound");
    soplex->changeUpperReal(i, value != INF ? value : soplex::infinity);
  }

  SoplexLp::Value SoplexLp::_getColUpperBound(int i) const {
    double value = soplex->upperReal(i);
    return value != soplex::infinity ? value : INF;
  }

  void SoplexLp::_setRowLowerBound(int i, Value lb) {
    LEMON_ASSERT(lb != INF, "Invalid bound");
    soplex->changeRangeReal(i, lb != -INF ? lb : -soplex::infinity, soplex->rhsReal(i));
  }

  SoplexLp::Value SoplexLp::_getRowLowerBound(int i) const {
    double res = soplex->lhsReal(i);
    return res == -soplex::infinity ? -INF : res;
  }

  void SoplexLp::_setRowUpperBound(int i, Value ub) {
    LEMON_ASSERT(ub != -INF, "Invalid bound");
    soplex->changeRangeReal(i, soplex->lhsReal(i), ub != INF ? ub : soplex::infinity);
  }

  SoplexLp::Value SoplexLp::_getRowUpperBound(int i) const {
    double res = soplex->rhsReal(i);
    return res == soplex::infinity ? INF : res;
  }

  void SoplexLp::_setObjCoeffs(ExprIterator b, ExprIterator e) {
    for (int j = 0; j < soplex->numColsReal(); ++j) {
      soplex->changeObjReal(j, 0.0);
    }
    for (ExprIterator it = b; it != e; ++it) {
      soplex->changeObjReal(it->first, it->second);
    }
  }

  void SoplexLp::_getObjCoeffs(InsertIterator b) const {
    for (int j = 0; j < soplex->numColsReal(); ++j) {
      Value coef = soplex->objReal(j);
      if (coef != 0.0) {
        *b = std::make_pair(j, coef);
        ++b;
      }
    }
  }

  void SoplexLp::_setObjCoeff(int i, Value objReal_coef) {
    soplex->changeObjReal(i, objReal_coef);
  }

  SoplexLp::Value SoplexLp::_getObjCoeff(int i) const {
    return soplex->objReal(i);
  }

  SoplexLp::SolveExitStatus SoplexLp::_solve() {

    _clear_temporals();

    _applyMessageLevel();

    soplex::SPxSolver::Status status = soplex->solve();

    switch (status) {
    case soplex::SPxSolver::OPTIMAL:
    case soplex::SPxSolver::INFEASIBLE:
    case soplex::SPxSolver::UNBOUNDED:
      return SOLVED;
    default:
      return UNSOLVED;
    }
  }

  SoplexLp::Value SoplexLp::_getPrimal(int i) const {
    if (_primal_values.empty()) {
      _primal_values.resize(soplex->numColsReal());
      soplex::Vector pv(_primal_values.size(), &_primal_values.front());
      soplex->getPrimal(pv);
    }
    return _primal_values[i];
  }

  SoplexLp::Value SoplexLp::_getDual(int i) const {
    if (_dual_values.empty()) {
      _dual_values.resize(soplex->numRowsReal());
      soplex::Vector dv(_dual_values.size(), &_dual_values.front());
      soplex->getDual(dv);
    }
    return _dual_values[i];
  }

  SoplexLp::Value SoplexLp::_getPrimalValue() const {
    return soplex->objValueReal();
  }

  SoplexLp::VarStatus SoplexLp::_getColStatus(int i) const {
    switch (soplex->basisColStatus(i)) {
    case soplex::SPxSolver::BASIC:
      return BASIC;
    case soplex::SPxSolver::ON_UPPER:
      return UPPER;
    case soplex::SPxSolver::ON_LOWER:
      return LOWER;
    case soplex::SPxSolver::FIXED:
      return FIXED;
    case soplex::SPxSolver::ZERO:
      return FREE;
    default:
      LEMON_ASSERT(false, "Wrong column status");
      return VarStatus();
    }
  }

  SoplexLp::VarStatus SoplexLp::_getRowStatus(int i) const {
    switch (soplex->basisRowStatus(i)) {
    case soplex::SPxSolver::BASIC:
      return BASIC;
    case soplex::SPxSolver::ON_UPPER:
      return UPPER;
    case soplex::SPxSolver::ON_LOWER:
      return LOWER;
    case soplex::SPxSolver::FIXED:
      return FIXED;
    case soplex::SPxSolver::ZERO:
      return FREE;
    default:
      LEMON_ASSERT(false, "Wrong row status");
      return VarStatus();
    }
  }

  SoplexLp::Value SoplexLp::_getPrimalRay(int i) const {
    if (_primal_ray.empty()) {
      _primal_ray.resize(soplex->numColsReal());
      soplex::Vector pv(_primal_ray.size(), &_primal_ray.front());
      soplex->getDualFarkas(pv);
    }
    return _primal_ray[i];
  }

  SoplexLp::Value SoplexLp::_getDualRay(int i) const {
    if (_dual_ray.empty()) {
      _dual_ray.resize(soplex->numRowsReal());
      soplex::Vector dv(_dual_ray.size(), &_dual_ray.front());
      soplex->getDualFarkas(dv);
    }
    return _dual_ray[i];
  }

  SoplexLp::ProblemType SoplexLp::_getPrimalType() const {
    switch (soplex->status()) {
    case soplex::SPxSolver::OPTIMAL:
      return OPTIMAL;
    case soplex::SPxSolver::UNBOUNDED:
      return UNBOUNDED;
    case soplex::SPxSolver::INFEASIBLE:
      return INFEASIBLE;
    default:
      return UNDEFINED;
    }
  }

  SoplexLp::ProblemType SoplexLp::_getDualType() const {
    switch (soplex->status()) {
    case soplex::SPxSolver::OPTIMAL:
      return OPTIMAL;
    case soplex::SPxSolver::UNBOUNDED:
      return UNBOUNDED;
    case soplex::SPxSolver::INFEASIBLE:
      return INFEASIBLE;
    default:
      return UNDEFINED;
    }
  }

  void SoplexLp::_setSense(Sense sense) {
      // STUB
  }

  SoplexLp::Sense SoplexLp::_getSense() const {
      return SoplexLp::Sense();
  }

  void SoplexLp::_clear() {
    soplex->clearBasis();
    _col_names.clear();
    _col_names_ref.clear();
    _row_names.clear();
    _row_names_ref.clear();
    cols.clear();
    rows.clear();
    _clear_temporals();
  }

  void SoplexLp::_messageLevel(MessageLevel level) {
    switch (level) {
    case MESSAGE_NOTHING:
      _message_level = -1;
      break;
    case MESSAGE_ERROR:
      _message_level = soplex::SPxOut::Verbosity::VERB_ERROR;
      break;
    case MESSAGE_WARNING:
      _message_level = soplex::SPxOut::Verbosity::VERB_WARNING;
      break;
    case MESSAGE_NORMAL:
      _message_level = soplex::SPxOut::Verbosity::VERB_INFO2;
      break;
    case MESSAGE_VERBOSE:
      _message_level = soplex::SPxOut::Verbosity::VERB_DEBUG;
      break;
    }
  }

  void SoplexLp::_applyMessageLevel() {
    soplex->setIntParam(soplex::SoPlex::VERBOSITY, _message_level);
  }

} //namespace lemon
