//==============================================================================
//
//	Copyright (c) 2013-
//	Authors:
//	* Dave Parker <d.a.parker@cs.bham.ac.uk> (University of Birmingham/Oxford)
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

import common.FastUtils;
import de.tum.in.naturals.set.NatBitSet;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.IntIterable;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntIterators;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.ints.IntSets;
import parser.State;
import parser.Values;
import prism.Pair;
import prism.PrismException;
import prism.PrismNotSupportedException;
import strat.MDStrategy;

import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.function.IntToDoubleFunction;

/**
* Explicit-state representation of a DTMC, constructed (implicitly)
* from an MDP and a memoryless deterministic (MD) adversary.
* This class is read-only: most of the data is pointers to other model info.
*/
public class DTMCFromMDPAndMDStrategy extends DTMCExplicit
{
	// Parent MDP
	protected MDP mdp;
	// MD strategy
	protected MDStrategy strat;

	/**
	 * Constructor: create from MDP and memoryless adversary.
	 */
	public DTMCFromMDPAndMDStrategy(MDP mdp, MDStrategy strat)
	{
		this.mdp = mdp;
		this.numStates = mdp.getNumStates();
		this.strat = strat;
	}

	@Override
	public void buildFromPrismExplicit(String filename) throws PrismException
	{
		throw new PrismNotSupportedException("Not supported");
	}

	// Accessors (for Model)

	@Override public IntIterable getInitialStates()
	{
		return mdp.getInitialStates();
	}

	@Override public boolean isInitialState(int i)
	{
		return mdp.isInitialState(i);
	}

	@Override public boolean isDeadlockState(int i)
	{
		return !strat.isChoiceDefined(i) || mdp.isDeadlockState(i);
	}

	@Override public List<State> getStatesList()
	{
		return mdp.getStatesList();
	}

	@Override public Values getConstantValues()
	{
		return mdp.getConstantValues();
	}

	@Override public int getNumTransitions()
	{
		int numTransitions = 0;
		for (int state = 0; state < numStates; state++) {
			if (strat.isChoiceDefined(state)) {
				numTransitions += mdp.getNumTransitions(state, strat.getChoiceIndex(state));
			}
		}
		return numTransitions;
	}

	@Override
	public Iterator<Int2ObjectMap.Entry<Pair<Double, Object>>> getTransitionsAndActionsIterator(int s)
	{
		if (strat.isChoiceDefined(s)) {
			int choiceIndex = strat.getChoiceIndex(s);
			return new DTMCExplicit.AddDefaultActionToTransitionsIterator(mdp.getTransitions(s, choiceIndex).iterator(), mdp.getAction(s, choiceIndex));
		} else {
			return Collections.emptyIterator();
		}
	}

	@Override public IntIterator getSuccessorsIterator(int s)
	{
		if (strat.isChoiceDefined(s)) {
			return mdp.getSuccessorsIterator(s, strat.getChoiceIndex(s));
		} else {
			return IntIterators.EMPTY_ITERATOR;
		}
	}

	@Override public boolean isSuccessor(int s1, int s2)
	{
		return strat.isChoiceDefined(s1) && mdp.someSuccessorsInSet(s1, strat.getChoiceIndex(s1), IntSets.singleton(s2));
	}

	@Override public boolean allSuccessorsInSet(int s, IntSet set)
	{
		return mdp.allSuccessorsInSet(s, strat.getChoiceIndex(s), set);
	}

	@Override public boolean someSuccessorsInSet(int s, IntSet set)
	{
		return mdp.someSuccessorsInSet(s, strat.getChoiceIndex(s), set);
	}

	@Override public void findDeadlocks(boolean fix) throws PrismException
	{
		throw new UnsupportedOperationException("Immutable object");
	}

	// Accessors (for DTMC)

	@Override public int getNumTransitions(int s)
	{
		return strat.isChoiceDefined(s) ? mdp.getNumTransitions(s, strat.getChoiceIndex(s)) : 0;
	}

	@Override public Iterator<Int2DoubleMap.Entry> getTransitionsIterator(int s)
	{
		if (strat.isChoiceDefined(s)) {
			return mdp.getTransitions(s, strat.getChoiceIndex(s)).iterator();
		} else {
			return Collections.emptyIterator();
		}
	}

	@Override public void prob0step(IntSet subset, IntSet u, NatBitSet result)
	{
		FastUtils.forEach(subset, numStates, i -> {
			if (!strat.isChoiceDefined(i)) {
				result.clear(i);
			} else {
				result.set(i, mdp.someSuccessorsInSet(i, strat.getChoiceIndex(i), u));
			}
		});
	}

	@Override public void prob1step(IntSet subset, IntSet u, IntSet v, NatBitSet result)
	{
		FastUtils.forEach(subset, numStates, i -> {
			if (!strat.isChoiceDefined(i)) {
				result.clear(i);
			} else {
				int j = strat.getChoiceIndex(i);
				result.set(i, mdp.someSuccessorsInSet(i, j, v) && mdp.allSuccessorsInSet(i, j, u));
			}
		});
	}

	@Override
	public double mvMultSingle(int s, IntToDoubleFunction vect)
	{
		return strat.isChoiceDefined(s) ? mdp.mvMultSingle(s, strat.getChoiceIndex(s), vect) : 0;
	}

	@Override
	public double mvMultJacSingle(int s, double vect[])
	{
		return strat.isChoiceDefined(s) ? mdp.mvMultJacSingle(s, strat.getChoiceIndex(s), vect) : 0;
	}

	@Override
	public boolean equals(Object o)
	{
		throw new RuntimeException("Not implemented yet");
	}
}
