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

import de.tum.in.naturals.set.NatBitSet;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.IntIterable;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntSet;
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
 * from an MDP and a memoryless adversary, specified as an array of integer indices.
 * This class is read-only: most of data is pointers to other model info.
 */
public class DTMCFromMDPMemorylessAdversary extends DTMCExplicit
{
	// Parent MDP
	protected MDP mdp;
	// Adversary (array of choice indices; -1 denotes no choice)
	protected MDStrategy adv;

	/**
	 * Constructor: create from MDP and memoryless adversary.
	 */
	public DTMCFromMDPMemorylessAdversary(MDP mdp, MDStrategy adv)
	{
		this.mdp = mdp;
		this.numStates = mdp.getNumStates();
		this.adv = adv;
	}

	@Override
	public void buildFromPrismExplicit(String filename) throws PrismException
	{
		throw new PrismNotSupportedException("Not supported");
	}

	// Accessors (for Model)

	@Override public int getNumStates()
	{
		return mdp.getNumStates();
	}

	@Override public int getNumInitialStates()
	{
		return mdp.getNumInitialStates();
	}

	@Override public IntIterable getInitialStates()
	{
		return mdp.getInitialStates();
	}

	@Override public int getFirstInitialState()
	{
		return mdp.getFirstInitialState();
	}

	@Override public boolean isInitialState(int i)
	{
		return mdp.isInitialState(i);
	}

	@Override public boolean isDeadlockState(int i)
	{
		return mdp.isDeadlockState(i);
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
		for (int s = 0; s < numStates; s++)
			if (adv.isChoiceDefined(s))
				numTransitions += mdp.getNumTransitions(s, adv.getChoiceIndex(s));
		return numTransitions;
	}

	@Override public IntIterator getSuccessorsIterator(final int s)
	{
		int choiceId = adv.getChoiceIndex(s);
		if (choiceId < 0)
			throw new RuntimeException("Choice not set for the queried state");
		else
			return mdp.getSuccessorsIterator(s, adv.getChoiceIndex(s));
	}

	@Override public boolean isSuccessor(int s1, int s2)
	{
		throw new RuntimeException("Not implemented yet");
	}

	@Override public boolean allSuccessorsInSet(int s, IntSet set)
	{
		int choiceId = adv.getChoiceIndex(s);
		if (choiceId < 0) {
			throw new RuntimeException("Choice not set for the queried state");
		}
		return mdp.getSuccessorsIntStream(s, choiceId).allMatch(set::contains);
	}

	@Override public boolean someSuccessorsInSet(int s, IntSet set)
	{
		int choiceId = adv.getChoiceIndex(s);
		if (choiceId < 0) {
			throw new RuntimeException("Choice not set for the queried state");
		}
		return mdp.getSuccessorsIntStream(s, choiceId).anyMatch(set::contains);
	}

	public int getNumChoices(int s)
	{
		// Always 1 for a DTMC
		return 1;
	}

	@Override public void findDeadlocks(boolean fix) throws PrismException
	{
		// No deadlocks by definition
	}

	@Override public void checkForDeadlocks() throws PrismException
	{
		// No deadlocks by definition
	}

	@Override public void checkForDeadlocks(IntSet except) throws PrismException
	{
		// No deadlocks by definition
	}

	@Override
	public String infoString()
	{
		return mdp.infoString() + " + " + "???"; // TODO
	}

	@Override
	public String infoStringTable()
	{
		return mdp.infoString() + " + " + "???\n"; // TODO
	}

	// Accessors (for DTMC)

	@Override
	public int getNumTransitions(int s)
	{
		return adv.isChoiceDefined(s) ? mdp.getNumTransitions(s, adv.getChoiceIndex(s)) : 0;
	}

	@Override
	public Iterator<Int2DoubleMap.Entry> getTransitionsIterator(int s)
	{
		if (adv.isChoiceDefined(s)) {
			return mdp.getTransitions(s, adv.getChoiceIndex(s)).iterator();
		} else {
			// Empty iterator
			return Collections.emptyIterator();
		}
	}

	@Override
	public Iterator<Int2ObjectMap.Entry<Pair<Double, Object>>> getTransitionsAndActionsIterator(int s)
	{
		if (adv.isChoiceDefined(s)) {
			return new DTMCExplicit.AddDefaultActionToTransitionsIterator(mdp.getTransitions(s, adv.getChoiceIndex(s)).iterator(), mdp.getAction(s, adv.getChoiceIndex(s)));
		} else {
			// Empty iterator
			return Collections.emptyIterator();
		}
	}

	@Override public void prob0step(IntSet subset, IntSet u, NatBitSet result)
	{
		// TODO
		throw new Error("Not yet supported");
	}

	@Override public void prob1step(IntSet subset, IntSet u, IntSet v, NatBitSet result)
	{
		for (int s = 0; s < numStates; s++) {
			if (subset.contains(s)) {
				result.set(s, mdp.prob1stepSingle(s, adv.getChoiceIndex(s), u, v));
			}
		}
	}

	@Override
	public double mvMultSingle(int s, IntToDoubleFunction vect)
	{
		return adv.isChoiceDefined(s) ? mdp.mvMultSingle(s, adv.getChoiceIndex(s), vect) : 0;
	}

	@Override
	public double mvMultJacSingle(int s, double vect[])
	{
		return adv.isChoiceDefined(s) ? mdp.mvMultJacSingle(s, adv.getChoiceIndex(s), vect) : 0;
	}

	@Override
	public boolean equals(Object o)
	{
		throw new RuntimeException("Not implemented yet");
	}
}
