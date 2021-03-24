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
import explicit.rewards.MCRewards;
import it.unimi.dsi.fastutil.ints.AbstractInt2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleAVLTreeMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.IntIterable;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntIterators;
import it.unimi.dsi.fastutil.ints.IntSet;
import parser.State;
import parser.Values;
import prism.PrismException;
import prism.PrismNotSupportedException;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.function.IntToDoubleFunction;

/**
 * Simple explicit-state representation of a DTMC, constructed (implicitly) as the embedded DTMC of a CTMC.
 * i.e. P(i,j) = R(i,j) / E(i) if E(i) > 0 and P(i,i) = 1 otherwise
 * where E(i) is the exit rate for state i: sum_j R(i,j).
 * This class is read-only: most of data is pointers to other model info.
 */
public class DTMCEmbeddedSimple extends DTMCExplicit
{
	// Parent CTMC
	protected CTMCSimple ctmc;
	// Exit rates vector
	protected double exitRates[];
	// Number of extra transitions added (just for stats)
	protected int numExtraTransitions;

	/**
	 * Constructor: create from CTMC.
	 */
	public DTMCEmbeddedSimple(CTMCSimple ctmc)
	{
		this.ctmc = ctmc;
		this.numStates = ctmc.getNumStates();
		// TODO: should we copy other stuff across too?
		exitRates = new double[numStates];
		numExtraTransitions = 0;
		for (int i = 0; i < numStates; i++) {
			exitRates[i] = ctmc.getTransitions(i).sum();
			if (exitRates[i] == 0)
				numExtraTransitions++;
		}
	}

	@Override
	public void buildFromPrismExplicit(String filename) throws PrismException
	{
		throw new PrismNotSupportedException("Not supported");
	}

	// Accessors (for Model)

	@Override public int getNumStates()
	{
		return ctmc.getNumStates();
	}

	@Override public int getNumInitialStates()
	{
		return ctmc.getNumInitialStates();
	}

	@Override public IntIterable getInitialStates()
	{
		return ctmc.getInitialStates();
	}

	@Override public int getFirstInitialState()
	{
		return ctmc.getFirstInitialState();
	}

	@Override public boolean isInitialState(int i)
	{
		return ctmc.isInitialState(i);
	}

	@Override public boolean isDeadlockState(int i)
	{
		return ctmc.isDeadlockState(i);
	}

	@Override public List<State> getStatesList()
	{
		return ctmc.getStatesList();
	}

	@Override public Values getConstantValues()
	{
		return ctmc.getConstantValues();
	}

	@Override public int getNumTransitions()
	{
		return ctmc.getNumTransitions() + numExtraTransitions;
	}

	@Override public IntIterator getSuccessorsIterator(int s)
	{
		if (exitRates[s] == 0) {
			return IntIterators.singleton(s);
		} else {
			return ctmc.getSuccessorsIterator(s);
		}
	}

	@Override public boolean isSuccessor(int s1, int s2)
	{
		return exitRates[s1] == 0 ? (s1 == s2) : ctmc.isSuccessor(s1, s2);
	}

	@Override public boolean allSuccessorsInSet(int s, IntSet set)
	{
		return exitRates[s] == 0 ? set.contains(s) : ctmc.allSuccessorsInSet(s, set);
	}

	@Override public boolean someSuccessorsInSet(int s, IntSet set)
	{
		return exitRates[s] == 0 ? set.contains(s) : ctmc.someSuccessorsInSet(s, set);
	}

	@Override public NatBitSet getLabelStates(String name)
	{
		return ctmc.getLabelStates(name);
	}

	@Override public boolean hasLabel(String name)
	{
		return ctmc.hasLabel(name);
	}

	@Override public Set<String> getLabels()
	{
		return ctmc.getLabels();
	}

	@Override public void addLabel(String name, NatBitSet states)
	{
		throw new RuntimeException("Can not add label to DTMCEmbeddedSimple");
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
		String s = "";
		s += numStates + " states (" + getNumInitialStates() + " initial)";
		s += ", " + getNumTransitions() + " transitions (incl. " + numExtraTransitions + " self-loops)";
		return s;
	}

	@Override
	public String infoStringTable()
	{
		String s = "";
		s += "States:      " + numStates + " (" + getNumInitialStates() + " initial)\n";
		s += "Transitions: " + getNumTransitions() + "\n";
		return s;
	}

	// Accessors (for DTMC)

	@Override public int getNumTransitions(int s)
	{
		if (exitRates[s] == 0) {
			return 1;
		} else {
			return ctmc.getNumTransitions(s);
		}
	}

	@Override public Iterator<Int2DoubleMap.Entry> getTransitionsIterator(int s)
	{
		if (exitRates[s] == 0) {
			// return prob-1 self-loop
			Int2DoubleMap m = new Int2DoubleAVLTreeMap();
			m.put(s, 1.0);
			return m.int2DoubleEntrySet().iterator();
		} else {
			Iterator<Int2DoubleMap.Entry> ctmcIterator = ctmc.getTransitionsIterator(s);

			// return iterator over entries, with probabilities divided by exitRates[s]
			double er = exitRates[s];
			return new Iterator<Int2DoubleMap.Entry>() {
				@Override
				public boolean hasNext()
				{
					return ctmcIterator.hasNext();
				}

				@Override
				public Int2DoubleMap.Entry next()
				{
					Int2DoubleMap.Entry ctmcEntry = ctmcIterator.next();
					return new AbstractInt2DoubleMap.BasicEntry(ctmcEntry.getIntKey(), ctmcEntry.getDoubleValue() / er);
				}

				@Override
				public void remove()
				{
					throw new UnsupportedOperationException();
				}
			};
		}
	}

	@Override public double mvMultSingle(int s, IntToDoubleFunction vect) {
		int k;
		double d, er, prob;
		Distribution distr;

		distr = ctmc.getTransitions(s);
		d = 0.0;
		er = exitRates[s];
		// Exit rate 0: prob 1 self-loop
		if (er == 0) {
			d += vect.applyAsDouble(s);
		}
		// Exit rate > 0
		else {
			d = distr.sumWeighted(vect);
			d /= er;
		}

		return d;
	}

	@Override
	public double mvMultJacSingle(int s, double vect[])
	{
		int k;
		double diag, d, er, prob;
		Distribution distr;

		distr = ctmc.getTransitions(s);
		diag = d = 0.0;
		er = exitRates[s];
		// Exit rate 0: prob 1 self-loop
		if (er == 0) {
			return 0.0;
		}
		// Exit rate > 0
		else {
			// (sum_{j!=s} P(s,j)*vect[j]) / (1-P(s,s))
			// = (sum_{j!=s} (R(s,j)/E(s))*vect[j]) / (1-(P(s,s)/E(s)))
			// = (sum_{j!=s} R(s,j)*vect[j]) / (E(s)-P(s,s))
			for (Int2DoubleMap.Entry e : distr) {
				k = e.getIntKey();
				prob = e.getDoubleValue();
				// Non-diagonal entries only
				if (k != s) {
					d += prob * vect[k];
				} else {
					diag = prob;
				}
			}
			d /= (er - diag);
		}

		return d;
	}

	//@Override
	public double mvMultRewJacSingle(int s, double vect[], MCRewards mcRewards)
	{
		int k;
		double diag, d, er, prob;
		Distribution distr;

		distr = ctmc.getTransitions(s);
		diag = d = 0.0;
		er = exitRates[s];
		// Exit rate 0: prob 1 self-loop
		if (er == 0) {
			return mcRewards.getStateReward(s);
		}
		// Exit rate > 0
		// (rew(s) + sum_{j!=s} P(s,j)*vect[j]) / (1-P(s,s))
		// = (rew(s) + sum_{j!=s} (R(s,j)/E(s))*vect[j]) / (1-(P(s,s)/E(s)))
		// = (E(s)*rew(s) + sum_{j!=s} R(s,j)*vect[j]) / (E(s)-P(s,s))
		d = er * mcRewards.getStateReward(s);
		for (Int2DoubleMap.Entry e : distr) {
			k = e.getIntKey();
			prob = e.getDoubleValue();
			// Non-diagonal entries only
			if (k != s) {
				d += prob * vect[k];
			} else {
				diag = prob;
			}
		}
		d /= (er - diag);

		return d;
	}

	@Override
	public void vmMult(double vect[], double result[])
	{
		int i, j;
		double prob, er;
		Distribution distr;

		// Initialise result to 0
		for (j = 0; j < numStates; j++) {
			result[j] = 0;
		}
		// Go through matrix elements (by row)
		for (i = 0; i < numStates; i++) {
			distr = ctmc.getTransitions(i);
			er = exitRates[i];
			// Exit rate 0: prob 1 self-loop
			if (er == 0) {
				result[i] += vect[i];
			}
			// Exit rate > 0
			else {
				for (Int2DoubleMap.Entry e : distr) {
					j = e.getIntKey();
					prob = e.getDoubleValue();
					result[j] += (prob / er) * vect[i];
				}
			}
		}
	}

	@Override
	public String toString()
	{
		int i, numStates;
		boolean first;
		StringBuilder s = new StringBuilder();
		s.append("ctmc: ").append(ctmc);
		first = true;
		s.append(", exitRates: [ ");
		numStates = getNumStates();
		for (i = 0; i < numStates; i++) {
			if (first)
				first = false;
			else
				s.append(", ");
			s.append(i).append(": ").append(exitRates[i]);
		}
		s.append(" ]");
		return s.toString();
	}

	@Override
	public boolean equals(Object o)
	{
		if (o == null || !(o instanceof DTMCEmbeddedSimple))
			return false;
		DTMCEmbeddedSimple dtmc = (DTMCEmbeddedSimple) o;
		if (!ctmc.equals(dtmc.ctmc))
			return false;
		if (!Arrays.equals(exitRates, dtmc.exitRates))
			return false;
		if (numExtraTransitions != dtmc.numExtraTransitions)
			return false;
		return true;
	}
}
