//==============================================================================
//
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
//	* Christian von Essen <christian.vonessen@imag.fr> (Verimag, Grenoble)
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
import common.IntDoubleConsumer;
import de.tum.in.naturals.set.NatBitSet;
import explicit.rewards.MDPRewards;
import it.unimi.dsi.fastutil.ints.Int2DoubleAVLTreeMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import parser.State;
import prism.PrismException;
import prism.PrismUtils;
import strat.MDStrategy;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.function.IntConsumer;
import java.util.function.IntToDoubleFunction;

/**
 * Sparse matrix (non-mutable) explicit-state representation of an MDP.
 * This is much faster to access than e.g. MDPSimple and should also be more compact.
 * The catch is that you have to create the model all in one go and then can't modify it.
 */
public class MDPSparse extends MDPExplicit
{
	// Sparse matrix storing transition function (Steps)
	/** Probabilities for each transition (array of size numTransitions) */
	protected double nonZeros[];
	/** Column (destination) indices for each transition (array of size numTransitions) */
	protected int cols[];
	/** Indices into nonZeros/cols giving the start of the transitions for each choice (distribution);
	 * array is of size numDistrs+1 and last entry is always equal to numTransitions */
	protected int choiceStarts[];
	/** Indices into choiceStarts giving the start of the choices for each state;
	 * array is of size numStates+1 and last entry is always equal to numDistrs */
	protected int rowStarts[];

	// Action labels
	/** Array of action labels for choices;
	 * if null, there are no actions; otherwise, is an array of size numDistrs */
	protected Object actions[];

	// Other statistics
	protected int numDistrs;
	protected int numTransitions;
	protected int maxNumDistrs;

	// Constructors

	/**
	 * Copy constructor (from MDPSimple).
	 */
	public MDPSparse(MDPSimple mdp) {
		initialise(mdp.getNumStates());
		copyFrom(mdp);
		// Copy stats
		numDistrs = mdp.getNumChoices();
		numTransitions = mdp.getNumTransitions();
		maxNumDistrs = mdp.getMaxNumChoices();

		nonZeros = new double[numTransitions];
		cols = new int[numTransitions];
		choiceStarts = new int[numDistrs + 1];
		rowStarts = new int[numStates + 1];
		actions = mdp.actions == null ? null : new Object[numDistrs];
		int j = 0;
		int k = 0;
		for (int i = 0; i < numStates; i++) {
			rowStarts[i] = j;
			if (mdp.actions != null) {
				int n = mdp.getNumChoices(i);
				for (int l = 0; l < n; l++) {
					actions[j + l] = mdp.getAction(i, l);
				}
			}

			for (Distribution distr : mdp.trans.get(i)) {
				choiceStarts[j] = k;
				for (Int2DoubleMap.Entry e : distr) {
					cols[k] = e.getIntKey();
					nonZeros[k] = e.getDoubleValue();
					k++;
				}
				// k += distr.keys.length;
				j++;
			}
		}

		choiceStarts[numDistrs] = numTransitions;
		rowStarts[numStates] = numDistrs;
	}

	/**
	 * Copy constructor (from MDPSimple). Optionally, transitions within choices
	 * are sorted (by ascending order of column index). Also, optionally, a state
	 * index permutation can be provided, i.e. old state index i becomes index permut[i].
	 * Note: a states list, if present, will not be permuted and should be set
	 * separately afterwards if required.
	 * @param mdp The MDP to copy
	 * @param sort Whether or not to sort column indices
	 * @param permut State space permutation
	 */
	public MDPSparse(MDPSimple mdp, boolean sort, int permut[])
	{
		int i, j, k, n;
		Int2DoubleMap sorted = null;
		int permutInv[];
		initialise(mdp.getNumStates());
		copyFrom(mdp, permut);
		// Copy stats
		numDistrs = mdp.getNumChoices();
		numTransitions = mdp.getNumTransitions();
		maxNumDistrs = mdp.getMaxNumChoices();
		// Compute the inverse of the permutation
		permutInv = new int[numStates];
		for (i = 0; i < numStates; i++) {
			permutInv[permut[i]] = i;
		}
		// Copy transition function
		if (sort) {
			sorted = new Int2DoubleAVLTreeMap();
		}
		nonZeros = new double[numTransitions];
		cols = new int[numTransitions];
		choiceStarts = new int[numDistrs + 1];
		rowStarts = new int[numStates + 1];
		actions = mdp.actions == null ? null : new Object[numDistrs];
		j = 0;
		k = 0;
		for (i = 0; i < numStates; i++) {
			rowStarts[i] = j;
			if (mdp.actions != null) {
				n = mdp.getNumChoices(permutInv[i]);
				for (int l = 0; l < n; l++) {
					actions[j + l] = mdp.getAction(permutInv[i], l);
				}
			}
			for (Distribution distr : mdp.trans.get(permutInv[i])) {
				choiceStarts[j] = k;
				for (Int2DoubleMap.Entry e : distr) {
					if (sort) {
						sorted.put(permut[e.getIntKey()], e.getDoubleValue());
					} else {
						cols[k] = permut[e.getIntKey()];
						nonZeros[k] = e.getDoubleValue();
						k++;
					}
				}
				if (sort) {
					for (Int2DoubleMap.Entry e : sorted.int2DoubleEntrySet()) {
						cols[k] = e.getIntKey();
						nonZeros[k] = e.getDoubleValue();
						k++;
					}
					sorted.clear();
				}
				j++;
			}
		}
		choiceStarts[numDistrs] = numTransitions;
		rowStarts[numStates] = numDistrs;
	}

	/**
	 * Copy constructor for a (sub-)MDP from a given MDP.
	 * The states and actions will be indexed as given by the order
	 * of the lists {@code states} and {@code actions}.
	 * @param mdp MDP to copy from
	 * @param states States to copy
	 * @param actions Actions to copy
	 */
	public MDPSparse(MDP mdp, IntList states, List<IntList> actions)
	{
		initialise(states.size());
		mdp.getInitialStates().forEach((IntConsumer) this::addInitialState);
		mdp.getDeadlockStates().forEach((IntConsumer) this::addDeadlockState);
		statesList = new ArrayList<>(states.size());
		List<State> statesList = mdp.getStatesList();
		states.forEach((IntConsumer) s -> statesList.add(statesList.get(s)));
		numDistrs = 0;
		numTransitions = 0;
		maxNumDistrs = 0;
		states.forEach((IntConsumer) state -> {
			IntList stateActions = actions.get(state);
			int numChoices = stateActions.size();
			numDistrs += numChoices;
			if (numChoices > maxNumDistrs) {
				maxNumDistrs = numChoices;
			}
			IntIterator iterator1 = stateActions.iterator();
			while (iterator1.hasNext()) {
				numTransitions += mdp.getNumTransitions(state, iterator1.nextInt());
			}
		});
		nonZeros = new double[numTransitions];
		cols = new int[numTransitions];
		choiceStarts = new int[numDistrs + 1];
		rowStarts = new int[numStates + 1];
		this.actions = new Object[numDistrs];
		int choiceIndex = 0;
		int colIndex = 0;
		int[] reverseStates = new int[mdp.getNumStates()];
		for (int i = 0; i < states.size(); i++) {
			reverseStates[states.getInt(i)] = i;
		}
		for (int i = 0; i < states.size(); i++) {
			int s = states.getInt(i);
			rowStarts[i] = choiceIndex;
			IntIterator iterator = actions.get(s).iterator();
			while (iterator.hasNext()) {
				int a = iterator.nextInt();
				choiceStarts[choiceIndex] = colIndex;
				this.actions[choiceIndex] = mdp.getAction(s, a);
				choiceIndex++;
				for (Int2DoubleMap.Entry next : mdp.getTransitions(s, a)) {
					cols[colIndex] = reverseStates[next.getIntKey()];
					nonZeros[colIndex] = next.getDoubleValue();
					colIndex++;
				}
			}
		}
		choiceStarts[numDistrs] = numTransitions;
		rowStarts[numStates] = numDistrs;
	}

	// Mutators (other)

	@Override
	public void initialise(int numStates)
	{
		super.initialise(numStates);
		numDistrs = 0;
		numTransitions = 0;
		maxNumDistrs = 0;
		actions = null;
	}

	@Override
	public void buildFromPrismExplicit(String filename) throws PrismException
	{
		String s, ss[];
		int i, j, k, iLast, kLast, jCount, kCount, n, lineNum = 0;
		double prob;

		// Open file
		try (BufferedReader in = new BufferedReader(new FileReader(new File(filename)))) {
			// Parse first line to get num states
			s = in.readLine();
			lineNum = 1;
			if (s == null) {
				in.close();
				throw new PrismException("Missing first line of .tra file");
			}
			ss = s.split(" ");
			n = Integer.parseInt(ss[0]);
			// Initialise
			initialise(n);
			// Set initial state (assume 0)
			initialStates.add(0);
			// Store stats
			numDistrs = Integer.parseInt(ss[1]);
			numTransitions = Integer.parseInt(ss[2]);
			// Go though list of transitions in file
			iLast = -1;
			kLast = -1;
			jCount = 0;
			kCount = 0;
			s = in.readLine();
			lineNum++;
			while (s != null) {
				s = s.trim();
				if (!s.isEmpty()) {
					ss = s.split(" ");
					i = Integer.parseInt(ss[0]);
					k = Integer.parseInt(ss[1]);
					j = Integer.parseInt(ss[2]);
					prob = Double.parseDouble(ss[3]);
					// For a new state
					if (i != iLast) {
						rowStarts[i] = kCount;
					}
					// For a new state or distribution
					if (i != iLast || k != kLast) {
						choiceStarts[kCount] = jCount;
						kCount++;
					}
					// Store transition
					cols[jCount] = j;
					nonZeros[jCount] = prob;
					// Prepare for next iter
					iLast = i;
					kLast = k;
					jCount++;
				}
				s = in.readLine();
				lineNum++;
			}
			choiceStarts[numDistrs] = numTransitions;
			rowStarts[numStates] = numDistrs;
			// Compute maxNumDistrs
			maxNumDistrs = 0;
			for (i = 0; i < numStates; i++) {
				maxNumDistrs = Math.max(maxNumDistrs, getNumChoices(i));
			}
			// Close file
			in.close();
			// Sanity checks
			if (kCount != numDistrs) {
				throw new PrismException("Choice count is wrong in tra file (" + kCount + "!=" + numTransitions + ")");
			}
			if (jCount != numTransitions) {
				throw new PrismException("Transition count is wrong in tra file (" + kCount + "!=" + numTransitions + ")");
			}
		} catch (IOException e) {
			System.exit(1);
		} catch (NumberFormatException e) {
			throw new PrismException("Problem in .tra file (line " + lineNum + ") for " + getModelType());
		}
	}

	// Accessors (for Model)

	@Override
	public int getNumTransitions()
	{
		return numTransitions;
	}

	@Override
	public IntIterator getSuccessorsIterator(final int s)
	{
		// Need to build set to avoid duplicates
		// So not necessarily the fastest method to access successors
		int start = choiceStarts[rowStarts[s]];
		int end = choiceStarts[rowStarts[s + 1]];
		IntSet succs = new IntOpenHashSet();
		for (int i = start; i < end; i++) {
			succs.add(cols[i]);
		}

		return succs.iterator();
	}

	@Override
	public boolean isSuccessor(int s1, int s2)
	{
		int j, k, l1, h1, l2, h2;
		l1 = rowStarts[s1];
		h1 = rowStarts[s1 + 1];
		for (j = l1; j < h1; j++) {
			l2 = choiceStarts[j];
			h2 = choiceStarts[j + 1];
			for (k = l2; k < h2; k++) {
				// Assume that only non-zero entries are stored
				if (cols[k] == s2) {
					return true;
				}
			}
		}
		return false;
	}

	@Override
	public boolean allSuccessorsInSet(int s, IntSet set)
	{
		int j, k, l1, h1, l2, h2;
		l1 = rowStarts[s];
		h1 = rowStarts[s + 1];
		for (j = l1; j < h1; j++) {
			l2 = choiceStarts[j];
			h2 = choiceStarts[j + 1];
			for (k = l2; k < h2; k++) {
				// Assume that only non-zero entries are stored
				if (!set.contains(cols[k])) {
					return false;
				}
			}
		}
		return true;
	}

	@Override
	public boolean someSuccessorsInSet(int s, IntSet set)
	{
		int j, k, l1, h1, l2, h2;
		l1 = rowStarts[s];
		h1 = rowStarts[s + 1];
		for (j = l1; j < h1; j++) {
			l2 = choiceStarts[j];
			h2 = choiceStarts[j + 1];
			for (k = l2; k < h2; k++) {
				// Assume that only non-zero entries are stored
				if (set.contains(cols[k])) {
					return true;
				}
			}
		}
		return false;
	}

	@Override
	public void findDeadlocks(boolean fix) throws PrismException
	{
		for (int i = 0; i < numStates; i++) {
			// Note that no distributions is a deadlock, not an empty distribution
			if (getNumChoices(i) == 0) {
				addDeadlockState(i);
				if (fix) {
					throw new PrismException("Can't fix deadlocks in an MDPSparse since it cannot be modified after construction");
				}
			}
		}
	}

	@Override
	public void checkForDeadlocks(IntSet except) throws PrismException
	{
		for (int i = 0; i < numStates; i++) {
			if (getNumChoices(i) == 0 && (except == null || !except.contains(i)))
				throw new PrismException("MDP has a deadlock in state " + i);
		}
	}

	// Accessors (for NondetModel)

	@Override
	public int getNumChoices(int s)
	{
		return rowStarts[s + 1] - rowStarts[s];
	}

	@Override
	public int getMaxNumChoices()
	{
		return maxNumDistrs;
	}

	@Override
	public int getNumChoices()
	{
		return numDistrs;
	}

	@Override
	public Object getAction(int s, int i)
	{
		return i < 0 || actions == null ? null : actions[rowStarts[s] + i];
	}

	@Override
	public boolean allSuccessorsInSet(int s, int i, IntSet set)
	{
		int j, k, l2, h2;
		j = rowStarts[s] + i;
		l2 = choiceStarts[j];
		h2 = choiceStarts[j + 1];
		for (k = l2; k < h2; k++) {
			// Assume that only non-zero entries are stored
			if (!set.contains(cols[k])) {
				return false;
			}
		}
		return true;
	}

	@Override
	public boolean someSuccessorsInSet(int s, int i, IntSet set)
	{
		int j, k, l2, h2;
		j = rowStarts[s] + i;
		l2 = choiceStarts[j];
		h2 = choiceStarts[j + 1];
		for (k = l2; k < h2; k++) {
			// Assume that only non-zero entries are stored
			if (set.contains(cols[k])) {
				return true;
			}
		}
		return false;
	}

	@Override
	public IntIterator getSuccessorsIterator(final int s, final int i)
	{
		int start = choiceStarts[rowStarts[s] + i];
		int end = choiceStarts[rowStarts[s] + i + 1];
		IntList successors = new IntArrayList(end - start);
		for (; start < end; start++) {
			successors.add(cols[start]);
		}
		return successors.iterator();
	}

	// Accessors (for MDP)

	@Override
	public int getNumTransitions(int s, int i)
	{
		return choiceStarts[rowStarts[s] + i + 1] - choiceStarts[rowStarts[s] + i];
	}

	class SparseEntry implements Int2DoubleMap.Entry {
		int i;

		@Override
		@Deprecated
		public Integer getKey()
		{
			return getIntKey();
		}

		@Override
		public int getIntKey() {
			return cols[i];
		}

		@Override
		@Deprecated
		public Double getValue()
		{
			return getDoubleValue();
		}

		@Override
		public double setValue(double value) {
			throw new UnsupportedOperationException();
		}

		@Override
		public double getDoubleValue() {
			return nonZeros[i];
		}

		@Deprecated
		@Override
		public Double setValue(Double value)
		{
			return setValue(value.doubleValue());
		}
	}

	class SparseIterator implements Iterator<Int2DoubleMap.Entry> {
		final int end;
		final SparseEntry entry;

		SparseIterator(int s, int i) {
			entry = new SparseEntry();
			entry.i = choiceStarts[rowStarts[s] + i] - 1;
			end = choiceStarts[rowStarts[s] + i + 1];
		}

		@Override
		public boolean hasNext()
		{
			return entry.i + 1 < end;
		}

		@Override
		public Int2DoubleMap.Entry next() {
			entry.i++;
			return entry;
		}

		@Override
		public void remove()
		{
			throw new UnsupportedOperationException();
		}
	}

	@Override
	public Iterable<Int2DoubleMap.Entry> getTransitions(final int s, final int i)
	{
		return () -> new SparseIterator(s, i);
	}

	@Override
	public void forEachTransition(int s, int i, IntDoubleConsumer consumer) {
		int current = choiceStarts[rowStarts[s] + i];
		int end = choiceStarts[rowStarts[s] + i + 1];

		for (; current < end; current++) {
			consumer.accept(cols[current], nonZeros[current]);
		}
	}

	@Override
	public void forEachTransition(int s, int i, IntConsumer consumer) {
		int current = choiceStarts[rowStarts[s] + i];
		int end = choiceStarts[rowStarts[s] + i + 1];

		for (; current < end; current++) {
			consumer.accept(cols[current]);
		}
	}

	@Override
	public void prob0step(IntSet subset, IntSet u, boolean forall, NatBitSet result)
	{
		int j, k, l1, h1, l2, h2;
		boolean b1, some;
		IntIterator iterator = subset.iterator();
		while (iterator.hasNext()) {
			int i = iterator.nextInt();
			b1 = forall; // there exists or for all
			l1 = rowStarts[i];
			h1 = rowStarts[i + 1];
			for (j = l1; j < h1; j++) {
				some = false;
				l2 = choiceStarts[j];
				h2 = choiceStarts[j + 1];
				for (k = l2; k < h2; k++) {
					// Assume that only non-zero entries are stored
					if (u.contains(cols[k])) {
						some = true;
						break;
					}
				}
				if (forall) {
					if (!some) {
						b1 = false;
						break;
					}
				} else {
					if (some) {
						b1 = true;
						break;
					}
				}
			}
			result.set(i, b1);
		}
	}

	@Override
	public void prob1Astep(IntSet subset, IntSet u, IntSet v, NatBitSet result)
	{
		int j, k, l1, h1, l2, h2;
		boolean b1, some, all;
		IntIterator iterator = subset.iterator();
		while (iterator.hasNext()) {
			int i = iterator.nextInt();
			b1 = true;
			l1 = rowStarts[i];
			h1 = rowStarts[i + 1];
			for (j = l1; j < h1; j++) {
				some = false;
				all = true;
				l2 = choiceStarts[j];
				h2 = choiceStarts[j + 1];
				for (k = l2; k < h2; k++) {
					// Assume that only non-zero entries are stored
					if (!u.contains(cols[k])) {
						all = false;
						break; // Stop early (already know b1 will be set to false)
					}

					if (!some) {
						some = v.contains(cols[k]);
					}
				}
				if (!(some && all)) {
					b1 = false;
					break;
				}
			}
			result.set(i, b1);
		}
	}

	@Override
	public void prob1Estep(IntSet subset, IntSet u, IntSet v, NatBitSet result, MDStrategy strat)
	{
		int l1, h1, l2, h2, stratCh = -1;
		boolean b1, some, all;
		IntIterator iterator = subset.iterator();
		while (iterator.hasNext()) {
			int i = iterator.nextInt();
			b1 = false;
			l1 = rowStarts[i];
			h1 = rowStarts[i + 1];
			for (int j = l1; j < h1; j++) {
				some = false;
				all = true;
				l2 = choiceStarts[j];
				h2 = choiceStarts[j + 1];
				for (int k = l2; k < h2; k++) {
					// Assume that only non-zero entries are stored
					if (!u.contains(cols[k])) {
						all = false;
						break; // Stop early (already know b1 will not be set to true)
					}

					if (!some) {
						some = v.contains(cols[k]);
					}
				}
				if (some && all) {
					b1 = true;
					// If strategy generation is enabled, remember optimal choice
					if (strat != null)
						stratCh = j - l1;
					break;
				}
			}
			// If strategy generation is enabled, store optimal choice
			// (only if this the first time we add the state to S^yes)
			if (strat != null & b1 & !result.contains(i)) {
				strat.setChoiceIndex(i, stratCh);
			}
			// Store result
			result.set(i, b1);
		}
	}

	@Override
	public void prob1step(IntSet subset, IntSet u, IntSet v, boolean forall, NatBitSet result)
	{
		IntIterator iterator = subset.iterator();
		while (iterator.hasNext()) {
			int i = iterator.nextInt();
			boolean b1 = forall; // there exists or for all
			int l1 = rowStarts[i];
			int h1 = rowStarts[i + 1];
			for (int j = l1; j < h1; j++) {
				boolean some = false;
				boolean all = true;
				int l2 = choiceStarts[j];
				int h2 = choiceStarts[j + 1];

				for (int k = l2; k < h2; k++) {
					// Assume that only non-zero entries are stored
					if (!some) {
						some = v.contains(cols[k]);
					}

					if (all) {
						all = u.contains(cols[k]);
					}
				}

				if (forall) {
					if (!(some && all)) {
						b1 = false;
						break;
					}
				} else {
					if (some && all) {
						b1 = true;
						break;
					}
				}
			}

			result.set(i, b1);
		}
	}

	@Override
	public boolean prob1stepSingle(int s, int i, IntSet u, IntSet v)
	{
		int j, k, l2, h2;
		boolean some, all;

		j = rowStarts[s] + i;
		some = false;
		all = true;
		l2 = choiceStarts[j];
		h2 = choiceStarts[j + 1];
		for (k = l2; k < h2; k++) {
			// Assume that only non-zero entries are stored
			if (v.contains(cols[k])) {
				some = true;
			}
			if (!u.contains(cols[k])) {
				all = false;
			}
		}
		return some && all;
	}

	@Override
	public double mvMultMinMaxSingle(int s, double vect[], boolean min, MDStrategy strat)
	{
		int j, k, l1, h1, l2, h2, stratCh = -1;
		double d, minmax;
		boolean first;

		minmax = 0;
		first = true;
		l1 = rowStarts[s];
		h1 = rowStarts[s + 1];
		for (j = l1; j < h1; j++) {
			// Compute sum for this distribution
			d = 0.0;
			l2 = choiceStarts[j];
			h2 = choiceStarts[j + 1];
			for (k = l2; k < h2; k++) {
				d += nonZeros[k] * vect[cols[k]];
			}
			// Check whether we have exceeded min/max so far
			if (first || (min && d < minmax) || (!min && d > minmax)) {
				minmax = d;
				// If strategy generation is enabled, remember optimal choice
				if (strat != null)
					stratCh = j - l1;
			}
			first = false;
		}
		// If strategy generation is enabled, store optimal choice
		if (strat != null & !first) {
			// For max, only remember strictly better choices
			if (min) {
				strat.setChoiceIndex(s, stratCh);
			} else if (!strat.isChoiceDefined(s) || minmax > vect[s]) {
				strat.setChoiceIndex(s, stratCh);
			}
		}

		return minmax;
	}

	@Override
	public IntList mvMultMinMaxSingleChoices(int s, double vect[], boolean min, double val)
	{
		int j, k, l1, h1, l2, h2;
		double d;

		// Create data structures to store strategy
		IntList res = new IntArrayList();
		// One row of matrix-vector operation
		l1 = rowStarts[s];
		h1 = rowStarts[s + 1];
		for (j = l1; j < h1; j++) {
			// Compute sum for this distribution
			d = 0.0;
			l2 = choiceStarts[j];
			h2 = choiceStarts[j + 1];
			for (k = l2; k < h2; k++) {
				d += nonZeros[k] * vect[cols[k]];
			}
			// Store strategy info if value matches
			if (PrismUtils.doublesAreClose(val, d, 1e-12, false)) {
				res.add(j - l1);
			}
		}

		return res;
	}

	@Override
	public double mvMultSingle(int s, int i, IntToDoubleFunction vect)
	{
		int j, k, l2, h2;
		double d;

		j = rowStarts[s] + i;
		// Compute sum for this distribution
		d = 0.0;
		l2 = choiceStarts[j];
		h2 = choiceStarts[j + 1];
		for (k = l2; k < h2; k++) {
			d += nonZeros[k] * vect.applyAsDouble(cols[k]);
		}

		return d;
	}

	@Override
	public double mvMultJacMinMaxSingle(int s, double vect[], boolean min, MDStrategy strat)
	{
		int j, k, l1, h1, l2, h2, stratCh = -1;
		double diag, d, minmax;
		boolean first;

		minmax = 0;
		first = true;
		l1 = rowStarts[s];
		h1 = rowStarts[s + 1];
		for (j = l1; j < h1; j++) {
			diag = 1.0;
			// Compute sum for this distribution
			d = 0.0;
			l2 = choiceStarts[j];
			h2 = choiceStarts[j + 1];
			for (k = l2; k < h2; k++) {
				if (cols[k] != s) {
					d += nonZeros[k] * vect[cols[k]];
				} else {
					diag -= nonZeros[k];
				}
			}
			if (diag > 0)
				d /= diag;
			// Check whether we have exceeded min/max so far
			if (first || (min && d < minmax) || (!min && d > minmax)) {
				minmax = d;
				// If strategy generation is enabled, remember optimal choice
				if (strat != null)
					stratCh = j - l1;
			}
			first = false;
		}
		// If strategy generation is enabled, store optimal choice
		if (strat != null & !first) {
			// For max, only remember strictly better choices
			if (min) {
				strat.setChoiceIndex(s, stratCh);
			} else if (!strat.isChoiceDefined(s) || minmax > vect[s]) {
				strat.setChoiceIndex(s, stratCh);
			}
		}

		return minmax;
	}

	@Override
	public double mvMultJacSingle(int s, int i, double vect[])
	{
		int j, k, l2, h2;
		double diag, d;

		j = rowStarts[s] + i;
		diag = 1.0;
		// Compute sum for this distribution
		d = 0.0;
		l2 = choiceStarts[j];
		h2 = choiceStarts[j + 1];
		for (k = l2; k < h2; k++) {
			if (cols[k] != s) {
				d += nonZeros[k] * vect[cols[k]];
			} else {
				diag -= nonZeros[k];
			}
		}
		if (diag > 0)
			d /= diag;

		return d;
	}

	@Override
	public double mvMultRewMinMaxSingle(int state, double vect[], MDPRewards mdpRewards, boolean min, MDStrategy strat)
	{
		int j, k, l1, h1, l2, h2, stratCh = -1;
		double d, minmax;
		boolean first;

		minmax = 0;
		first = true;
		l1 = rowStarts[state];
		h1 = rowStarts[state + 1];
		for (j = l1; j < h1; j++) {
			// Compute sum for this distribution
			d = mdpRewards.getTransitionReward(state, j - l1);
			l2 = choiceStarts[j];
			h2 = choiceStarts[j + 1];
			for (k = l2; k < h2; k++) {
				d += nonZeros[k] * vect[cols[k]];
			}
			// Check whether we have exceeded min/max so far
			if (first || (min && d < minmax) || (!min && d > minmax)) {
				minmax = d;
				// If strategy generation is enabled, remember optimal choice
				if (strat != null)
					stratCh = j - l1;
			}
			first = false;
		}
		// Add state reward (doesn't affect min/max)
		minmax += mdpRewards.getStateReward(state);
		// If strategy generation is enabled, store optimal choice
		if (strat != null & !first) {
			// For max, only remember strictly better choices
			if (min) {
				strat.setChoiceIndex(state, stratCh);
			} else if (!strat.isChoiceDefined(state) || minmax > vect[state]) {
				strat.setChoiceIndex(state, stratCh);
			}
		}

		return minmax;
	}

	/**
	 * See Puterman Aperiodicity Transformation (Section 8.5.4)
	 * v_i+1(s) = tau * (r(s) + r(s, a)) + (1 - tau) * v_i(s) + sum_{j \in S\{s}} ( tau * p(j | s, a) * v_i(j) )
	 * @param state Row index
	 * @param vect Vector to multiply by
	 * @param mdpRewards The rewards
	 * @param min Min or max for (true=min, false=max)
	 * @param strat Storage for (memoryless) strategy choice indices (ignored if null)
	 * @param tau Aperiodic transformation parameter
	 */
	@Override
	public double mvMultRewMinMaxAperiodicSingle(int state, IntSet allowedChoices, IntToDoubleFunction vect, MDPRewards mdpRewards,
			boolean min, MDStrategy strat, double tau)
	{
		assert tau > 0d;

		int strategyChoice = -1;
		double value = 0d;
		boolean first = true;

		final int choiceIndexStart = rowStarts[state];
		final int choiceIndexEnd = rowStarts[state + 1];
		final int numberOfChoices = choiceIndexEnd - choiceIndexStart;
		assert numberOfChoices == getNumChoices(state);

		IntIterator actionIterator = FastUtils.iterator(allowedChoices, numberOfChoices);
		while (actionIterator.hasNext()) {
			int currentActionIndex = actionIterator.nextInt();

			final int currentChoiceIndex = currentActionIndex + choiceIndexStart;
			final int successorDistributionStart = choiceStarts[currentChoiceIndex];
			final int successorDistributionEnd = choiceStarts[currentChoiceIndex + 1];

			// Compute sum for this distribution
			double choiceRewardSum = 0d;
			for (int distributionIndex = successorDistributionStart; distributionIndex < successorDistributionEnd; distributionIndex++) {
				final int successorState = cols[distributionIndex];
				final double transitionProbability = nonZeros[distributionIndex];
				choiceRewardSum += transitionProbability * vect.applyAsDouble(successorState);
			}
			choiceRewardSum *= tau;
			choiceRewardSum += ((1d - tau) * vect.applyAsDouble(state))
					+ tau * mdpRewards.getTransitionReward(state, currentChoiceIndex - choiceIndexStart);

			// Check whether we have exceeded min/max so far
			if (first || (min && choiceRewardSum < value) || (!min && choiceRewardSum > value)) {
				value = choiceRewardSum;
				// If strategy generation is enabled, remember optimal choice
				if (strat != null) {
					strategyChoice = currentChoiceIndex - choiceIndexStart;
				}
			}
			first = false;
		}

		// Add state reward (doesn't affect min/max)
		value += tau * mdpRewards.getStateReward(state);

		// If strategy generation is enabled, store optimal choice
		if (strat != null & !first) {
			// For max, only remember strictly better choices
			if (min) {
				strat.setChoiceIndex(state, strategyChoice);
			} else if (!strat.isChoiceDefined(state) || value > vect.applyAsDouble(state)) {
				strat.setChoiceIndex(state, strategyChoice);
			}
		}

		return value;
	}

	@Override
	public double mvMultRewJacMinMaxSingle(int s, double vect[], MDPRewards mdpRewards, boolean min, MDStrategy strat)
	{
		int j, k, l1, h1, l2, h2, stratCh = -1;
		double diag, d, minmax;
		boolean first;

		minmax = 0;
		first = true;
		l1 = rowStarts[s];
		h1 = rowStarts[s + 1];
		for (j = l1; j < h1; j++) {
			diag = 1.0;
			// Compute sum for this distribution
			// (note: have to add state rewards in the loop for Jacobi)
			d = mdpRewards.getStateReward(s);
			d += mdpRewards.getTransitionReward(s, j - l1);
			l2 = choiceStarts[j];
			h2 = choiceStarts[j + 1];
			for (k = l2; k < h2; k++) {
				if (cols[k] != s) {
					d += nonZeros[k] * vect[cols[k]];
				} else {
					diag -= nonZeros[k];
				}
			}
			if (diag > 0)
				d /= diag;
			// Catch special case of probability 1 self-loop (Jacobi does it wrong)
			if (h2 - l2 == 1 && cols[l2] == s) {
				d = Double.POSITIVE_INFINITY;
			}
			// Check whether we have exceeded min/max so far
			if (first || (min && d < minmax) || (!min && d > minmax)) {
				minmax = d;
				// If strategy generation is enabled, remember optimal choice
				if (strat != null)
					stratCh = j - l1;
			}
			first = false;
		}
		// If strategy generation is enabled, store optimal choice
		if (strat != null & !first) {
			// For max, only remember strictly better choices
			if (min) {
				strat.setChoiceIndex(s, stratCh);
			} else if (!strat.isChoiceDefined(s) || minmax > vect[s]) {
				strat.setChoiceIndex(s, stratCh);
			}
		}

		return minmax;
	}

	@Override
	public IntList mvMultRewMinMaxSingleChoices(int state, double vect[], MDPRewards mdpRewards, boolean min, double val)
	{
		int j, k, l1, h1, l2, h2;
		double d;

		// Create data structures to store strategy
		IntList res = new IntArrayList();
		// One row of matrix-vector operation
		l1 = rowStarts[state];
		h1 = rowStarts[state + 1];
		for (j = l1; j < h1; j++) {
			// Compute sum for this distribution
			d = mdpRewards.getTransitionReward(state, j - l1);
			l2 = choiceStarts[j];
			h2 = choiceStarts[j + 1];
			for (k = l2; k < h2; k++) {
				d += nonZeros[k] * vect[cols[k]];
			}
			d += mdpRewards.getStateReward(state);
			// Store strategy info if value matches
			if (PrismUtils.doublesAreCloseRel(val, d, 1e-12)) {
				res.add(j - l1);
			}
		}

		return res;
	}

	@Override
	public void mvMultRight(int[] states, MDStrategy strat, double[] source, double[] dest)
	{
		for (int s : states) {
			int j = rowStarts[s] + strat.getChoiceIndex(s);
			int l2 = choiceStarts[j];
			int h2 = choiceStarts[j + 1];
			for (int k = l2; k < h2; k++) {
				dest[cols[k]] += nonZeros[k] * source[s];
			}
		}
	}

	// Standard methods

	@Override
	public String toString()
	{
		int i, j, k, l1, h1, l2, h2;
		Object o;
		StringBuilder builder = new StringBuilder();
		builder.append("[ ");
		for (i = 0; i < numStates; i++) {
			if (i > 0)
				builder.append(", ");
			builder.append(i).append(": [");
			l1 = rowStarts[i];
			h1 = rowStarts[i + 1];
			for (j = l1; j < h1; j++) {
				if (j > l1)
					builder.append(",");
				o = getAction(i, j - l1);
				if (o != null)
					builder.append(o).append(":");
				builder.append("{");
				l2 = choiceStarts[j];
				h2 = choiceStarts[j + 1];
				for (k = l2; k < h2; k++) {
					if (k > l2)
						builder.append(", ");
					builder.append(cols[k]).append("=").append(nonZeros[k]);
				}
				builder.append("}");
			}
			builder.append("]");
		}
		builder.append(" ]");

		return builder.toString();
	}

	@Override
	public boolean equals(Object o)
	{
		if (o == null || !(o instanceof MDPSparse))
			return false;
		MDPSparse mdp = (MDPSparse) o;
		if (numStates != mdp.numStates)
			return false;
		if (!initialStates.equals(mdp.initialStates))
			return false;
		if (!Utils.doubleArraysAreEqual(nonZeros, mdp.nonZeros))
			return false;
		if (!Utils.intArraysAreEqual(cols, mdp.cols))
			return false;
		if (!Utils.intArraysAreEqual(choiceStarts, mdp.choiceStarts))
			return false;
		if (!Utils.intArraysAreEqual(rowStarts, mdp.rowStarts))
			return false;
		// TODO: compare actions (complicated: null = null,null,null,...)
		return true;
	}
}
