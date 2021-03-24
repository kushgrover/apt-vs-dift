//==============================================================================
//
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
//
//------------------------------------------------------------------------------
//
//	This file is part of PRISM.
//
//	PRISM is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//
//	PRISM is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with PRISM; if not, write to the Free Software Foundation,
//	Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//==============================================================================

package explicit;

import com.google.common.collect.Iterables;
import com.google.common.collect.Iterators;
import de.tum.in.naturals.set.NatBitSet;
import it.unimi.dsi.fastutil.ints.IntIterable;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.ints.IntSets;
import parser.State;
import parser.Values;
import parser.VarList;
import prism.ModelType;
import prism.PrismException;
import prism.PrismLog;

import java.io.File;
import java.util.List;
import java.util.Set;

/**
 * Interface for (abstract) classes that provide (read-only) access to an explicit-state model.
 */
public interface Model
{
	// Accessors

	/**
	 * Get the type of this model.
	 */
	ModelType getModelType();

	/**
	 * Get the number of states.
	 */
	int getNumStates();

	/**
	 * Get the number of initial states.
	 */
	default int getNumInitialStates()
	{
		return Iterables.size(getInitialStates());
	}

	/**
	 * Get iterator over initial state list.
	 */
	IntIterable getInitialStates();

	/**
	 * Get the index of the first initial state
	 * (i.e. the one with the lowest index).
	 * Returns -1 if there are no initial states.
	 */
	default int getFirstInitialState() {
		IntIterator iterator = getInitialStates().iterator();
		if (iterator.hasNext()) {
			return iterator.nextInt();
		}
		return -1;
	}

	/**
	 * Check whether a state is an initial state.
	 */
	default boolean isInitialState(int i) {
		IntIterator iterator = getInitialStates().iterator();
		while (iterator.hasNext()) {
			if (i == iterator.nextInt()) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Get the number of states that are/were deadlocks.
	 * (Such states may have been fixed at build-time by adding self-loops)
	 */
	default int getNumDeadlockStates() {
		return Iterables.size(getDeadlockStates());
	}

	/**
	 * Get iterator over states that are/were deadlocks.
	 * (Such states may have been fixed at build-time by adding self-loops)
	 */
	IntIterable getDeadlockStates();

	/**
	 * Get list of states that are/were deadlocks.
	 * (Such states may have been fixed at build-time by adding self-loops)
	 */
	StateValues getDeadlockStatesList();

	/**
	 * Get the index of the first state that is/was a deadlock.
	 * (i.e. the one with the lowest index).
	 * Returns -1 if there are no initial states.
	 */
	default int getFirstDeadlockState() {
		IntIterator iterator = getDeadlockStates().iterator();
		if (iterator.hasNext()) {
			return iterator.nextInt();
		}
		return -1;
	}

	/**
	 * Check whether a state is/was deadlock.
	 * (Such states may have been fixed at build-time by adding self-loops)
	 */
	default boolean isDeadlockState(int i) {
		IntIterator iterator = getDeadlockStates().iterator();
		while (iterator.hasNext()) {
			if (i == iterator.nextInt()) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Get access to a list of states (optionally stored).
	 */
	List<State> getStatesList();

	/**
	 * Get access to the VarList (optionally stored)
	 */
	VarList getVarList();

	/**
	 * Get access to a list of constant values (optionally stored).
	 */
	Values getConstantValues();

	/**
	 * Get the states that satisfy a label in this model (optionally stored).
	 * Returns null if there is no label of this name.
	 */
	NatBitSet getLabelStates(String name);

	/**
	 * Get the labels that are (optionally) stored.
	 * Returns an empty set if there are no labels.
	 */
	Set<String> getLabels();

	/**
	 * Returns true if a label with the given name is attached to this model
	 */
	boolean hasLabel(String name);

	/**
	 * Get the total number of transitions in the model.
	 */
	default int getNumTransitions() {
		int numTransitions = 0;
		int numStates = getNumStates();
		for (int state = 0; state < numStates; state++) {
			numTransitions += Iterators.size(getSuccessorsIterator(state));
		}
		return numTransitions;
	}

	/**
	 * Get an iterator over the successors of state s.
	 */
	IntIterator getSuccessorsIterator(int s);

	/**
	 * Returns true if state s2 is a successor of state s1.
	 */
	default boolean isSuccessor(int s1, int s2) {
		IntIterator successorsIterator = getSuccessorsIterator(s1);
		while (successorsIterator.hasNext()) {
			if (successorsIterator.nextInt() == s2) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Check if all the successor states of a state are in a set.
	 *
	 * @param s   The state to check
	 * @param set The set to test for inclusion
	 */
	default boolean allSuccessorsInSet(int s, IntSet set) {
		IntIterator successorsIterator = getSuccessorsIterator(s);
		while (successorsIterator.hasNext()) {
			if (!set.contains(successorsIterator.nextInt())) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Check if any successor states of a state are in a set.
	 *
	 * @param s   The state to check
	 * @param set The set to test for inclusion
	 */
	default boolean someSuccessorsInSet(int s, IntSet set) {
		IntIterator successorsIterator = getSuccessorsIterator(s);
		while (successorsIterator.hasNext()) {
			if (set.contains(successorsIterator.nextInt())) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Find all deadlock states and store this information in the model.
	 * If requested (if fix=true) and if needed (i.e. for DTMCs/CTMCs),
	 * fix deadlocks by adding self-loops in these states.
	 * The set of deadlocks (before any possible fixing) can be obtained from {@link #getDeadlockStates()}.
	 *
	 * @throws PrismException if the model is unable to fix deadlocks because it is non-mutable.
	 */
	void findDeadlocks(boolean fix) throws PrismException;

	/**
	 * Checks for deadlocks and throws an exception if any exist.
	 */
	default void checkForDeadlocks() throws PrismException {
		checkForDeadlocks(IntSets.EMPTY_SET);
	}

	/**
	 * Checks for deadlocks and throws an exception if any exist.
	 * States in 'except' (If non-null) are excluded from the check.
	 */
	default void checkForDeadlocks(IntSet except) throws PrismException {
		for (int i = 0; i < getNumStates(); i++) {
			if (except.contains(i)) {
				continue;
			}
			if (!getSuccessorsIterator(i).hasNext()) {
				throw new PrismException("Deadlock in state " + i);
			}
		}
	}

	/**
	 * Export to explicit format readable by PRISM (i.e. a .tra file, etc.).
	 */
	void exportToPrismExplicit(String baseFilename) throws PrismException;

	/**
	 * Export transition matrix to explicit format readable by PRISM (i.e. a .tra file).
	 */
	void exportToPrismExplicitTra(String filename) throws PrismException;

	/**
	 * Export transition matrix to explicit format readable by PRISM (i.e. a .tra file).
	 */
	void exportToPrismExplicitTra(File file) throws PrismException;

	/**
	 * Export transition matrix to explicit format readable by PRISM (i.e. a .tra file).
	 */
	void exportToPrismExplicitTra(PrismLog log);

	/**
	 * Export to a dot file.
	 *
	 * @param filename Name of file to export to
	 */
	void exportToDotFile(String filename) throws PrismException;

	/**
	 * Export to a dot file, highlighting states in 'mark'.
	 *
	 * @param filename Name of file to export to
	 * @param mark     States to highlight (ignored if null)
	 */
	void exportToDotFile(String filename, IntSet mark) throws PrismException;

	/**
	 * Export to a dot file.
	 *
	 * @param out PrismLog to export to
	 */
	void exportToDotFile(PrismLog out);

	/**
	 * Export to a dot file, highlighting states in 'mark'.
	 *
	 * @param out  PrismLog to export to
	 * @param mark States to highlight (ignored if null)
	 */
	void exportToDotFile(PrismLog out, IntSet mark);

	/**
	 * Export to a dot file, highlighting states in 'mark'.
	 *
	 * @param out        PrismLog to export to
	 * @param mark       States to highlight (ignored if null)
	 * @param showStates Show state info on nodes?
	 */
	void exportToDotFile(PrismLog out, IntSet mark, boolean showStates);

	/**
	 * Export to a equivalent PRISM language model description.
	 */
	void exportToPrismLanguage(String filename) throws PrismException;

	/**
	 * Export states list.
	 */
	void exportStates(int exportType, VarList varList, PrismLog log) throws PrismException;

	/**
	 * Report info/stats about the model as a string.
	 */
	String infoString();

	/**
	 * Report info/stats about the model, tabulated, as a string.
	 */
	String infoStringTable();

	/**
	 * Has this model a stored PredecessorRelation?
	 */
	boolean hasStoredPredecessorRelation();

	/**
	 * If there is a PredecessorRelation stored for this model, return that.
	 * Otherwise, create one and return that. If {@code storeIfNew},
	 * store it for later use.
	 *
	 * @param parent     a PrismComponent (for obtaining the log)
	 * @param storeIfNew if the predecessor relation is newly created, store it
	 */
	PredecessorRelation getPredecessorRelation(prism.PrismComponent parent, boolean storeIfNew);

	/**
	 * Clear any stored predecessor relation, e.g., because the model was modified
	 */
	void clearPredecessorRelation();

}
