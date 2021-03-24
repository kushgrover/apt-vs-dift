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
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntIterators;
import it.unimi.dsi.fastutil.ints.IntSet;
import prism.PrismException;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.function.IntConsumer;
import java.util.function.IntToDoubleFunction;

/**
 * Simple explicit-state representation of a DTMC.
 */
public class DTMCSimple extends DTMCExplicit implements ModelSimple
{
	// Transition matrix (distribution list)
	protected List<Distribution> trans;

	// Other statistics
	protected int numTransitions;

	// Constructors

	/**
	 * Constructor: empty DTMC.
	 */
	public DTMCSimple()
	{
		initialise(0);
	}

	/**
	 * Constructor: new DTMC with fixed number of states.
	 */
	public DTMCSimple(int numStates)
	{
		initialise(numStates);
	}

	/**
	 * Copy constructor.
	 */
	public DTMCSimple(DTMCSimple dtmc)
	{
		this(dtmc.numStates);
		copyFrom(dtmc);
		for (int i = 0; i < numStates; i++) {
			trans.set(i, new Distribution(dtmc.trans.get(i)));
		}
		numTransitions = dtmc.numTransitions;
	}

	/**
	 * Construct a DTMC from an existing one and a state index permutation,
	 * i.e. in which state index i becomes index permut[i].
	 * Pointer to states list is NOT copied (since now wrong).
	 * Note: have to build new Distributions from scratch anyway to do this,
	 * so may as well provide this functionality as a constructor.
	 */
	public DTMCSimple(DTMCSimple dtmc, int permut[])
	{
		this(dtmc.numStates);
		copyFrom(dtmc, permut);
		for (int i = 0; i < numStates; i++) {
			trans.set(permut[i], new Distribution(dtmc.trans.get(i), permut));
		}
		numTransitions = dtmc.numTransitions;
	}

	// Mutators (for ModelSimple)

	@Override
	public void initialise(int numStates)
	{
		super.initialise(numStates);
		trans = new ArrayList<>(numStates);
		for (int i = 0; i < numStates; i++) {
			trans.add(new Distribution());
		}
	}

	@Override
	public void clearState(int i)
	{
		// Do nothing if state does not exist
		if (i >= numStates || i < 0) {
			return;
		}
		// Clear data structures and update stats
		numTransitions -= trans.get(i).size();
		trans.get(i).clear();
	}

	@Override
	public int addState()
	{
		addStates(1);
		return numStates - 1;
	}

	@Override
	public void addStates(int numToAdd)
	{
		for (int i = 0; i < numToAdd; i++) {
			trans.add(new Distribution());
			numStates++;
		}
	}

	@Override
	public void buildFromPrismExplicit(String filename) throws PrismException
	{
		String s, ss[];
		int i, j, n, lineNum = 0;
		double prob;

		// Open file for reading, automatic close when done
		try (BufferedReader in = new BufferedReader(new FileReader(new File(filename)))) {
			// Parse first line to get num states
			s = in.readLine();
			lineNum = 1;
			if (s == null) {
				throw new PrismException("Missing first line of .tra file");
			}
			ss = s.split(" ");
			n = Integer.parseInt(ss[0]);
			// Initialise
			initialise(n);
			// Go though list of transitions in file
			s = in.readLine();
			lineNum++;
			while (s != null) {
				s = s.trim();
				if (!s.isEmpty()) {
					ss = s.split(" ");
					i = Integer.parseInt(ss[0]);
					j = Integer.parseInt(ss[1]);
					prob = Double.parseDouble(ss[2]);
					setProbability(i, j, prob);
				}
				s = in.readLine();
				lineNum++;
			}
		} catch (IOException e) {
			throw new PrismException("File I/O error reading from \"" + filename + "\": " + e.getMessage());
		} catch (NumberFormatException e) {
			throw new PrismException("Problem in .tra file (line " + lineNum + ") for " + getModelType());
		}
	}

	// Mutators (other)

	/**
	 * Set the probability for a transition.
	 */
	public void setProbability(int i, int j, double prob)
	{
		Distribution distr = trans.get(i);
		if (distr.get(j) != 0.0) {
			numTransitions--;
		}
		if (prob != 0.0) {
			numTransitions++;
		}
		distr.set(j, prob);
	}

	/**
	 * Add to the probability for a transition.
	 */
	public void addToProbability(int i, int j, double prob)
	{
		if (!trans.get(i).add(j, prob)) {
			if (prob != 0.0) {
				numTransitions++;
			}
		}
	}

	// Accessors (for Model)

	@Override
	public int getNumTransitions()
	{
		return numTransitions;
	}

	@Override
	public IntIterator getSuccessorsIterator(int s)
	{
		return IntIterators.wrap(trans.get(s).getSupport());
	}

	@Override
	public boolean isSuccessor(int s1, int s2)
	{
		return trans.get(s1).contains(s2);
	}

	@Override
	public boolean allSuccessorsInSet(int s, IntSet set)
	{
		return (trans.get(s).isSubsetOf(set));
	}

	@Override
	public boolean someSuccessorsInSet(int s, IntSet set)
	{
		return (trans.get(s).containsOneOf(set));
	}

	@Override
	public void findDeadlocks(boolean fix) throws PrismException
	{
		for (int i = 0; i < numStates; i++) {
			if (trans.get(i).isEmpty()) {
				addDeadlockState(i);
				if (fix)
					setProbability(i, i, 1.0);
			}
		}
	}

	@Override
	public void checkForDeadlocks(IntSet except) throws PrismException
	{
		for (int i = 0; i < numStates; i++) {
			if (trans.get(i).isEmpty() && (except == null || !except.contains(i)))
				throw new PrismException("DTMC has a deadlock in state " + i);
		}
	}

	// Accessors (for DTMC)

	@Override
	public int getNumTransitions(int s)
	{
		return trans.get(s).size();
	}

	@Override
	public Iterator<Int2DoubleMap.Entry> getTransitionsIterator(int s)
	{
		return trans.get(s).iterator();
	}

	@Override
	public void prob0step(IntSet subset, IntSet u, NatBitSet result)
	{
		subset.forEach((IntConsumer) i -> result.set(i, trans.get(i).containsOneOf(u)));
	}

	@Override
	public void prob1step(IntSet subset, IntSet u, IntSet v, NatBitSet result)
	{
		subset.forEach((IntConsumer) i -> {
			Distribution distr = trans.get(i);
			result.set(i, distr.containsOneOf(v) && distr.isSubsetOf(u));
		});
	}

	@Override
	public double mvMultSingle(int s, double vect[])
	{
		return trans.get(s).sumWeighted(vect);
	}

	@Override
	public double mvMultSingle(int s, IntToDoubleFunction vect)
	{
		return trans.get(s).sumWeighted(vect);
	}

	@Override
	public double mvMultJacSingle(int s, double vect[])
	{
		Distribution distr = trans.get(s);
		double diag = 1.0;
		double d = 0.0d;

		for (Int2DoubleMap.Entry e : distr) {
			int k = e.getIntKey();
			double prob = e.getDoubleValue();
			if (k != s) {
				d += prob * vect[k];
			} else {
				diag -= prob;
			}
		}
		if (diag > 0) {
			d /= diag;
		}

		return d;
	}

	@Override
	public void vmMult(double vect[], double result[])
	{
		// Initialise result to 0
		Arrays.fill(result, 0, numStates, 0d);
		// Go through matrix elements (by row)
		for (int i = 0; i < numStates; i++) {
			Distribution distr = trans.get(i);
			for (Int2DoubleMap.Entry e : distr) {
				int j = e.getIntKey();
				double prob = e.getDoubleValue();
				result[j] += prob * vect[i];
			}
		}
	}

	// Accessors (other)

	/**
	 * Get the transitions (a distribution) for state s.
	 */
	public Distribution getTransitions(int s)
	{
		return trans.get(s);
	}

	// Standard methods

	@Override
	public boolean equals(Object o)
	{
		if (o == null || !(o instanceof DTMCSimple))
			return false;
		if (!super.equals(o))
			return false;
		DTMCSimple dtmc = (DTMCSimple) o;
		if (!trans.equals(dtmc.trans))
			return false;
		if (numTransitions != dtmc.numTransitions)
			return false;
		return true;
	}
}
