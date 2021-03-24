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

import common.FastUtils;
import common.IntDoubleConsumer;
import de.tum.in.naturals.set.NatBitSet;
import explicit.rewards.MDPRewards;
import it.unimi.dsi.fastutil.doubles.DoubleComparator;
import it.unimi.dsi.fastutil.doubles.DoubleComparators;
import it.unimi.dsi.fastutil.ints.Int2DoubleAVLTreeMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleSortedMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntIterators;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntSet;
import prism.ModelType;
import prism.PrismException;
import prism.PrismLog;
import prism.PrismUtils;
import strat.MDStrategy;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.function.IntConsumer;
import java.util.function.IntToDoubleFunction;

import static com.google.common.base.Preconditions.checkArgument;

/**
 * Interface for classes that provide (read) access to an explicit-state MDP.
 */
public interface MDP extends NondetModel
{
	@Override
	default ModelType getModelType()
	{
		return ModelType.MDP;
	}

	default Type getMDPType()
	{
		return Type.UNKNOWN;
	}

	default Set<Info> getMDPInformation()
	{
		return Collections.emptySet();
	}

	@Override
	default boolean areAllChoiceActionsUnique()
	{
		int numStates = getNumStates();
		Set<Object> sActions = new HashSet<>();
		for (int s = 0; s < numStates; s++) {
			int n = getNumChoices(s);
			if (n > 1) {
				sActions.clear();
				for (int i = 0; i < n; i++) {
					if (!sActions.add(getAction(s, i))) {
						return false;
					}
				}
			}
		}
		return true;
	}

	/**
	 * Get an iterator over the transitions from choice {@code i} of state {@code s}.
	 */
	Iterable<Int2DoubleMap.Entry> getTransitions(int s, int i);

	@Override default IntIterator getSuccessorsIterator(int s) {
		int numChoices = getNumChoices(s);
		if (numChoices == 0) {
			return IntIterators.EMPTY_ITERATOR;
		}
		IntIterator[] iterators = new IntIterator[numChoices];
		for (int i = 0; i < iterators.length; i++) {
			iterators[i] = getSuccessorsIterator(s, i);
		}
		return IntIterators.concat(iterators);
	}

	default Iterator<Int2DoubleMap.Entry> getTransitionsIterator(int stateIndex, int choiceIndex)
	{
		return getTransitions(stateIndex, choiceIndex).iterator();
	}

	default void forEachTransition(int s, int i, IntDoubleConsumer consumer)
	{
		getTransitions(s, i).forEach(entry -> consumer.accept(entry.getIntKey(), entry.getDoubleValue()));
	}

	default void forEachChoice(int state, IntConsumer consumer)
	{
		int numChoices = getNumChoices(state);
		for (int choice = 0; choice < numChoices; choice++) {
			consumer.accept(choice);
		}
	}

	/**
	 * Perform a single step of precomputation algorithm Prob0, i.e., for states i in {@code subset},
	 * set bit i of {@code result} iff, for all/some choices,
	 * there is a transition to a state in {@code u}.
	 * Quantification over choices is determined by {@code forall}.
	 *
	 * @param subset Only compute for these states
	 * @param u      Set of states {@code u}
	 * @param forall For-all or there-exists (true=for-all, false=there-exists)
	 * @param result Store results here
	 */
	default void prob0step(IntSet subset, IntSet u, boolean forall, NatBitSet result)
	{
		if (forall) {
			prob0Astep(subset, u, result);
		} else {
			prob0Estep(subset, u, result);
		}
	}

	default void prob0Astep(IntSet subset, IntSet u, NatBitSet result)
	{
		subset.forEach((IntConsumer) state -> {
			// Don't process the state if it was already solved
			if (result.contains(state)) {
				return;
			}
			for (int choice = 0; choice < getNumChoices(state); choice++) {
				if (!someSuccessorsInSet(state, choice, u)) {
					result.clear(state);
					return;
				}
			}
			result.set(state, state);
		});
	}

	default void prob0Estep(IntSet subset, IntSet u, NatBitSet result)
	{
		subset.forEach((IntConsumer) state -> {
			// Don't process the state if it was already solved
			if (result.contains(state)) {
				return;
			}
			for (int choice = 0; choice < getNumChoices(state); choice++) {
				if (someSuccessorsInSet(state, choice, u)) {
					result.set(state);
					return;
				}
			}
			result.clear(state);
		});
	}

	/**
	 * Perform a single step of precomputation algorithm Prob1A, i.e., for states i in {@code subset},
	 * set bit i of {@code result} iff, for all choices,
	 * there is a transition to a state in {@code v} and all transitions go to states in {@code u}.
	 *
	 * @param subset Only compute for these states
	 * @param u      Set of states {@code u}
	 * @param v      Set of states {@code v}
	 * @param result Store results here
	 */
	default void prob1Astep(IntSet subset, IntSet u, IntSet v, NatBitSet result)
	{
		subset.forEach((IntConsumer) state -> {
			// Don't process the state if it was already solved
			if (result.contains(state)) {
				return;
			}
			for (int choice = 0; choice < getNumChoices(state); choice++) {
				if (!prob1stepSingle(state, choice, u, v)) {
					result.clear(state);
					return;
				}
			}
			result.set(state);
		});
	}

	/**
	 * Perform a single step of precomputation algorithm Prob1E, i.e., for states i in {@code subset},
	 * set bit i of {@code result} iff, for some choice,
	 * there is a transition to a state in {@code v} and all transitions go to states in {@code u}.
	 * Optionally, store optimal (memoryless) strategy info for 1 states.
	 *
	 * @param subset Only compute for these states
	 * @param u      Set of states {@code u}
	 * @param v      Set of states {@code v}
	 * @param result Store results here
	 * @param strat  Storage for (memoryless) strategy choice indices (ignored if null)
	 */
	default void prob1Estep(IntSet subset, IntSet u, IntSet v, NatBitSet result, MDStrategy strat)
	{
		subset.forEach((IntConsumer) state -> {
			// Don't process the state if it was already solved
			if (result.contains(state)) {
				return;
			}
			int stateChoices = getNumChoices(state);
			for (int choice = 0; choice < stateChoices; choice++) {
				if (prob1stepSingle(state, choice, u, v)) {
					// If strategy generation is enabled, store optimal choice
					if (strat != null) {
						strat.setChoiceIndex(state, choice);
					}
					result.set(state);
					return;
				}
			}
			result.clear(state);
		});
	}

	/**
	 * Perform a single step of precomputation algorithm Prob1, i.e., for states i in {@code subset},
	 * set bit i of {@code result} iff, for all/some choices,
	 * there is a transition to a state in {@code v} and all transitions go to states in {@code u}.
	 * Quantification over choices is determined by {@code forall}.
	 *
	 * @param subset Only compute for these states
	 * @param u      Set of states {@code u}
	 * @param v      Set of states {@code v}
	 * @param forall For-all or there-exists (true=for-all, false=there-exists)
	 * @param result Store results here
	 */
	default void prob1step(IntSet subset, IntSet u, IntSet v, boolean forall, NatBitSet result)
	{
		if (forall) {
			prob1Astep(subset, u, v, result);
		} else {
			prob1Estep(subset, u, v, result, null);
		}
	}

	/**
	 * Perform a single step of precomputation algorithm Prob1 for a single state/choice,
	 * i.e., return whether there is a transition to a state in {@code v} and all transitions go to states in {@code u}.
	 *
	 * @param s State (row) index
	 * @param i Choice index
	 * @param u Set of states {@code u}
	 * @param v Set of states {@code v}
	 */
	default boolean prob1stepSingle(int s, int i, IntSet u, IntSet v)
	{
		IntIterator iterator = getSuccessorsIterator(s, i);
		boolean hasV = false;
		while (iterator.hasNext()) {
			int successor = iterator.nextInt();
			if (!u.contains(successor)) {
				return false;
			}
			hasV |= v.contains(successor);
		}
		return hasV;
	}

	/**
	 * Do a single row of matrix-vector multiplication for a specific choice.
	 *
	 * @param state  State (row) index
	 * @param choice Choice index
	 * @param vect   Vector to multiply by
	 */
	default double mvMultSingle(int state, int choice, IntToDoubleFunction vect)
	{
		Iterator<Int2DoubleMap.Entry> transitionsIterator = getTransitionsIterator(state, choice);
		double sum = 0.0d;
		while (transitionsIterator.hasNext()) {
			Int2DoubleMap.Entry transition = transitionsIterator.next();
			sum += transition.getDoubleValue() + vect.applyAsDouble(transition.getIntKey());
		}
		return sum;
	}

	/**
	 * Do a single row of matrix-vector multiplication for a specific choice.
	 *
	 * @param state  State (row) index
	 * @param choice Choice index
	 * @param vect   Vector to multiply by
	 */
	default double mvMultSingle(int state, int choice, double[] vect)
	{
		Iterator<Int2DoubleMap.Entry> transitionsIterator = getTransitionsIterator(state, choice);
		double sum = 0.0d;
		while (transitionsIterator.hasNext()) {
			Int2DoubleMap.Entry transition = transitionsIterator.next();
			sum += transition.getDoubleValue() * vect[transition.getIntKey()];
		}
		return sum;
	}

	default double mvMultSingleAperiodic(int s, int i, IntToDoubleFunction vect, double tau)
	{
		double successorValues = 0.0d;
		for (Int2DoubleMap.Entry entry : getTransitions(s, i)) {
			successorValues += entry.getDoubleValue() * vect.applyAsDouble(entry.getIntKey());
		}
		if (tau == 1d) {
			return successorValues;
		}
		return (1 - tau) * vect.applyAsDouble(s) + tau * successorValues;
	}

	/**
	 * Do a single row of Jacobi-style matrix-vector multiplication for a specific choice.
	 * i.e. return min/max_k { (sum_{j!=s} P_k(s,j)*vect[j]) / 1-P_k(s,s) }
	 *
	 * @param s    Row index
	 * @param i    Choice index
	 * @param vect Vector to multiply by
	 */
	default double mvMultJacSingle(int s, int i, double vect[])
	{
		double diag = 1.0;
		double sum = 0.0;
		Iterator<Int2DoubleMap.Entry> transitionsIterator = getTransitionsIterator(s, i);

		while (transitionsIterator.hasNext()) {
			Int2DoubleMap.Entry e = transitionsIterator.next();
			int k = e.getIntKey();
			double prob = e.getDoubleValue();
			if (k != s) {
				sum += prob * vect[k];
			} else {
				diag -= prob;
			}
		}
		if (diag > 0d) {
			sum /= diag;
		}

		return sum;
	}

	/**
	 * Do a single row of matrix-vector multiplication and sum of rewards for a specific choice.
	 * i.e. rew(s) + rew_k(s) + sum_j P_k(s,j)*vect[j]
	 *
	 * @param s          State (row) index
	 * @param i          Choice index
	 * @param vect       Vector to multiply by
	 * @param mdpRewards The rewards
	 */
	default double mvMultRewSingle(int s, int i, double vect[], MDPRewards mdpRewards)
	{
		double sum = mvMultSingle(s, i, vect);
		sum += mdpRewards.getTransitionReward(s, i);
		sum += mdpRewards.getStateReward(s);
		return sum;
	}

	/**
	 * Do a single row of matrix-vector multiplication and sum of rewards for a specific choice.
	 * i.e. rew(s) + rew_k(s) + sum_j P_k(s,j)*vect[j]
	 *
	 * @param s          State (row) index
	 * @param i          Choice index
	 * @param vect       Vector to multiply by
	 * @param mdpRewards The rewards
	 */
	default double mvMultRewSingle(int s, int i, IntToDoubleFunction vect, MDPRewards mdpRewards)
	{
		double sum = mvMultSingle(s, i, vect);
		sum += mdpRewards.getTransitionReward(s, i);
		sum += mdpRewards.getStateReward(s);
		return sum;
	}

	/**
	 * Do a matrix-vector multiplication followed by min/max, i.e. one step of value iteration,
	 * i.e. for all s: result[s] = min/max_k { sum_j P_k(s,j)*vect[j] }
	 * Optionally, store optimal (memoryless) strategy info.
	 *
	 * @param vect   Vector to multiply by
	 * @param min    Min or max for (true=min, false=max)
	 * @param result Vector to store result in
	 * @param strat  Storage for (memoryless) strategy choice indices (ignored if null)
	 */
	default void mvMultMinMax(double vect[], boolean min, double result[], IntSet states, MDStrategy strat)
	{
		states.forEach((IntConsumer) s -> result[s] = mvMultMinMaxSingle(s, vect, min, strat));
	}

	/**
	 * Do a Gauss-Seidel-style matrix-vector multiplication followed by min/max.
	 * i.e. for all s: vect[s] = min/max_k { (sum_{j!=s} P_k(s,j)*vect[j]) / 1-P_k(s,s) }
	 * and store new values directly in {@code vect} as computed.
	 * The maximum (absolute/relative) difference between old/new
	 * elements of {@code vect} is also returned.
	 * Optionally, store optimal (memoryless) strategy info.
	 *
	 * @param vect  Vector to multiply by (and store the result in)
	 * @param min   Min or max for (true=min, false=max)
	 * @param strat Storage for (memoryless) strategy choice indices (ignored if null)
	 * @return The maximum difference between old/new elements of {@code vect}
	 */
	default double mvMultGSMinMax(double vect[], boolean min, IntSet states, boolean absolute, MDStrategy strat)
	{
		double d, diff, maxDiff = 0.0;

		IntIterator iterator = states.iterator();
		while (iterator.hasNext()) {
			int state = iterator.nextInt();
			d = mvMultJacMinMaxSingle(state, vect, min, strat);
			diff = absolute ? (Math.abs(d - vect[state])) : (Math.abs(d - vect[state]) / d);
			maxDiff = diff > maxDiff ? diff : maxDiff;
			vect[state] = d;
		}

		return maxDiff;
	}

	/**
	 * Do a matrix-vector multiplication and sum of rewards followed by min/max, i.e. one step of value iteration.
	 * i.e. for all s: result[s] = min/max_k { rew(s) + rew_k(s) + sum_j P_k(s,j)*vect[j] }
	 * Optionally, store optimal (memoryless) strategy info.
	 *
	 * @param vect       Vector to multiply by
	 * @param mdpRewards The rewards
	 * @param min        Min or max for (true=min, false=max)
	 * @param result     Vector to store result in
	 * @param strat      Storage for (memoryless) strategy choice indices (ignored if null)
	 */
	default void mvMultRewMinMax(double vect[], MDPRewards mdpRewards, boolean min, double result[], IntSet states, MDStrategy strat)
	{
		states.forEach((IntConsumer) s -> result[s] = mvMultRewMinMaxSingle(s, vect, mdpRewards, min, strat));
	}

	/**
	 * Do a matrix-vector multiplication and sum of rewards followed by min/max, i.e. one step of value iteration.
	 * This function uses the transformation from Puterman, section 8.5.4 while performing a single step
	 * of value iteration. It can be thought of as having a self-loop on every state so that the periodicity
	 * becomes 1.
	 * i.e. for all s: v_i+1(s) = min/max_a { tau * (r(s) + r(s, a)) + (1 - tau) * v_i(s) + sum_{j \in S\{s}} ( tau * p(j | s, a) * v_i(j) ) }
	 * Optionally, store optimal (memoryless) strategy info.
	 *
	 * @param stateValues Vector to multiply by
	 * @param mdpRewards  The rewards
	 * @param min         Min or max for (true=min, false=max)
	 * @param result      Vector to store result in
	 * @param mec         Only do multiplication for these rows (ignored if null)
	 * @param strategy    Storage for (memoryless) strategy choice indices (ignored if null)
	 * @param tau         Aperiodic transformation parameter
	 */
	default void mvMultRewMinMaxAperiodic(IntToDoubleFunction stateValues, MDPRewards mdpRewards, boolean min, double result[], MEC mec,
			MDStrategy strategy, double tau)
	{
		IntIterator iterator = FastUtils.iterator(mec, getNumStates());
		while (iterator.hasNext()) {
			int state = iterator.nextInt();
			IntSet allowedActions = mec == null ? null : mec.actions.get(state);
			result[state] = mvMultRewMinMaxAperiodicSingle(state, allowedActions, stateValues, mdpRewards, min, strategy, tau);
		}
	}

	/**
	 * Do a Gauss-Seidel-style matrix-vector multiplication and sum of rewards followed by min/max.
	 * i.e. for all s: vect[s] = min/max_k { rew(s) + rew_k(s) + (sum_{j!=s} P_k(s,j)*vect[j]) / 1-P_k(s,s) }
	 * and store new values directly in {@code vect} as computed.
	 * The maximum (absolute/relative) difference between old/new
	 * elements of {@code vect} is also returned.
	 * Optionally, store optimal (memoryless) strategy info.
	 *
	 * @param vect       Vector to multiply by (and store the result in)
	 * @param mdpRewards The rewards
	 * @param min        Min or max for (true=min, false=max)
	 * @param absolute   If true, compute absolute, rather than relative, difference
	 * @param strat      Storage for (memoryless) strategy choice indices (ignored if null)
	 * @return The maximum difference between old/new elements of {@code vect}
	 */
	default double mvMultRewGSMinMax(double vect[], MDPRewards mdpRewards, boolean min, IntSet states, boolean absolute, MDStrategy strat)
	{
		double maxDiff = 0.0;
		IntIterator iterator = states.iterator();
		while (iterator.hasNext()) {
			int s = iterator.nextInt();
			double d = mvMultRewJacMinMaxSingle(s, vect, mdpRewards, min, strat);
			double diff = absolute ? (Math.abs(d - vect[s])) : (Math.abs(d - vect[s]) / d);
			maxDiff = diff > maxDiff ? diff : maxDiff;
			vect[s] = d;
		}
		// Use this code instead for backwards Gauss-Seidel
		/*for (s = numStates - 1; s >= 0; s--) {
			if (subset.get(s)) {
				d = mvMultRewJacMinMaxSingle(s, vect, mdpRewards, min);
				diff = absolute ? (Math.abs(d - vect[s])) : (Math.abs(d - vect[s]) / d);
				maxDiff = diff > maxDiff ? diff : maxDiff;
				vect[s] = d;
			}
		}*/
		return maxDiff;
	}

	default double mvMultMinMaxSingleGeneric(int state, IntIterator choices, StateChoiceValueFunction stateValue, boolean min,
			IntToDoubleFunction stratImproveVect, MDStrategy strat)
	{
		if (choices == null) {
			choices = IntIterators.fromTo(0, getNumChoices(state));
		}
		if (!choices.hasNext()) {
			return 0d;
		}

		double optimalValue = min ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
		int optimalChoice = -1;
		DoubleComparator comparator = min ? DoubleComparators.OPPOSITE_COMPARATOR : DoubleComparators.NATURAL_COMPARATOR;

		while (choices.hasNext()) {
			int choice = choices.nextInt();
			// Compute sum for this choice
			double value = stateValue.value(state, choice);

			// Check whether we have exceeded min/max so far
			if (comparator.compare(optimalValue, value) < 0) {
				optimalValue = value;
				optimalChoice = choice;
			}
		}

		// If strategy generation is enabled, store optimal choice
		if (strat != null) {
			// For max, only remember strictly better choices
			if (min) {
				strat.setChoiceIndex(state, optimalChoice);
			} else if (!strat.isChoiceDefined(state) || optimalValue > stratImproveVect.applyAsDouble(state)) {
				strat.setChoiceIndex(state, optimalChoice);
			}
		}
		return optimalValue;
	}

	/**
	 * Do a single row of matrix-vector multiplication followed by min/max,
	 * i.e. return min/max_k { sum_j P_k(s,j)*vect[j] }
	 * Optionally, store optimal (memoryless) strategy info.
	 *
	 * @param s     Row index
	 * @param vect  Vector to multiply by
	 * @param min   Min or max for (true=min, false=max)
	 * @param strat Storage for (memoryless) strategy choice indices (ignored if null)
	 */
	default double mvMultMinMaxSingle(int s, double vect[], boolean min, MDStrategy strat)
	{
		StateChoiceValueFunction stateValue = (state, choice) -> mvMultSingle(state, choice, vect);
		return mvMultMinMaxSingleGeneric(s, null, stateValue, min, i -> vect[i], strat);
	}

	/**
	 * Determine which choices result in min/max after a single row of matrix-vector multiplication.
	 *
	 * @param s    Row index
	 * @param vect Vector to multiply by
	 * @param min  Min or max (true=min, false=max)
	 * @param val  Min or max value to match
	 */
	default IntList mvMultMinMaxSingleChoices(int s, double vect[], boolean min, double val)
	{
		// Create data structures to store strategy
		IntList res = new IntArrayList();
		// One row of matrix-vector operation
		forEachChoice(s, choice -> {
			// Compute sum for this distribution
			double sum = mvMultSingle(s, choice, vect);
			// Store strategy info if value matches
			if (PrismUtils.doublesAreEqual(val, sum)) {
				res.add(choice);
			}
		});

		return res;
	}

	/**
	 * Do a single row of Jacobi-style matrix-vector multiplication followed by min/max.
	 * i.e. return min/max_k { (sum_{j!=s} P_k(s,j)*vect[j]) / 1-P_k(s,s) }
	 * Optionally, store optimal (memoryless) strategy info.
	 *
	 * @param s     Row index
	 * @param vect  Vector to multiply by
	 * @param min   Min or max for (true=min, false=max)
	 * @param strat Storage for (memoryless) strategy choice indices (ignored if null)
	 */
	default double mvMultJacMinMaxSingle(int s, double vect[], boolean min, MDStrategy strat)
	{
		StateChoiceValueFunction stateValue = (state, choice) -> mvMultJacSingle(state, choice, vect);
		return mvMultMinMaxSingleGeneric(s, null, stateValue, min, i -> vect[i], strat);
	}

	/**
	 * Do a single row of matrix-vector multiplication and sum of rewards followed by min/max.
	 * i.e. return min/max_k { rew(s) + rew_k(s) + sum_j P_k(s,j)*vect[j] }
	 * Optionally, store optimal (memoryless) strategy info.
	 *
	 * @param state      Row index
	 * @param vect       Vector to multiply by
	 * @param mdpRewards The rewards
	 * @param min        Min or max for (true=min, false=max)
	 * @param strat      Storage for (memoryless) strategy choice indices (ignored if null)
	 */
	default double mvMultRewMinMaxSingle(int state, double vect[], MDPRewards mdpRewards, boolean min, MDStrategy strat)
	{
		StateChoiceValueFunction stateValue = (s, i) -> mvMultRewSingle(s, i, vect, mdpRewards);
		return mvMultMinMaxSingleGeneric(state, null, stateValue, min, i -> vect[i], strat);
	}

	/**
	 * Do a single row of matrix-vector multiplication and sum of rewards followed by min/max.
	 * This function uses the transformation from Puterman, section 8.5.4 while performing a single step
	 * of value iteration. It can be thought of as having a self-loop on every state so that the periodicity
	 * becomes 1.
	 * i.e. return min/max_a { tau * (r(s) + r(s, a)) + (1 - tau) * v_i(s) + sum_{j \in S\{s}} ( tau * p(j | s, a) * v_i(j) ) }
	 * Optionally, store optimal (memoryless) strategy info.
	 *
	 * @param state      Row index
	 * @param vect       Vector to multiply by
	 * @param mdpRewards The rewards
	 * @param min        Min or max for (true=min, false=max)
	 * @param strat      Storage for (memoryless) strategy choice indices (ignored if null)
	 * @param tau        Aperiodic transformation parameter
	 */
	default double mvMultRewMinMaxAperiodicSingle(int state, IntSet allowedChoices, IntToDoubleFunction vect, MDPRewards mdpRewards,
			boolean min, MDStrategy strat, double tau)
	{
		checkArgument(0 < tau && tau <= 1d, "Aperiodicity parameter has to be in (0, 1]");

		StateChoiceValueFunction stateValue = (s, i) -> {
			// Compute sum for this choice
			double sum = mvMultRewSingle(s, i, vect, mdpRewards);
			if (tau < 1d) {
				// Re-weigh with aperodicity parameter
				sum = ((1d - tau) * vect.applyAsDouble(state)) + tau * sum;
			}
			return sum;
		};

		return mvMultMinMaxSingleGeneric(state, allowedChoices == null ? null : allowedChoices.iterator(), stateValue, min, vect, strat);
	}

	/**
	 * Do a single row of Jacobi-style matrix-vector multiplication and sum of rewards followed by min/max.
	 * i.e. return min/max_k { rew(s) + rew_k(s) + (sum_{j!=s} P_k(s,j)*vect[j]) / 1-P_k(s,s) }
	 * Optionally, store optimal (memoryless) strategy info.
	 *
	 * @param s          State (row) index
	 * @param vect       Vector to multiply by
	 * @param mdpRewards The rewards
	 * @param min        Min or max for (true=min, false=max)
	 * @param strat      Storage for (memoryless) strategy choice indices (ignored if null)
	 */
	default double mvMultRewJacMinMaxSingle(int s, double vect[], MDPRewards mdpRewards, boolean min, MDStrategy strat)
	{
		StateChoiceValueFunction stateValue = (state, choice) -> {
			// (note: have to add state rewards in the loop for Jacobi)
			double diag = 1.0;
			double sum = mdpRewards.getStateReward(state) + mdpRewards.getTransitionReward(state, choice);
			Iterator<Int2DoubleMap.Entry> transitionsIterator = getTransitionsIterator(state, choice);

			boolean onlySelfLoop = true;
			while (transitionsIterator.hasNext()) {
				Int2DoubleMap.Entry e = transitionsIterator.next();
				int k = e.getIntKey();
				double prob = e.getDoubleValue();
				if (k == s) {
					diag -= prob;
				} else {
					onlySelfLoop = false;
					sum += prob * vect[k];
				}
			}
			if (onlySelfLoop) {
				sum = Double.POSITIVE_INFINITY;
			} else if (diag > 0d) {
				sum /= diag;
			}

			return sum;
		};
		return mvMultMinMaxSingleGeneric(s, null, stateValue, min, i -> vect[i], strat);
	}

	/**
	 * Determine which choices result in min/max after a single row of matrix-vector multiplication and sum of rewards.
	 *
	 * @param state      State (row) index
	 * @param vect       Vector to multiply by
	 * @param mdpRewards The rewards
	 * @param min        Min or max (true=min, false=max)
	 * @param val        Min or max value to match
	 */
	default IntList mvMultRewMinMaxSingleChoices(int state, double vect[], MDPRewards mdpRewards, boolean min, double val)
	{
		// Create data structures to store choices
		IntList res = new IntArrayList();
		// One row of matrix-vector operation

		for (int choice = 0; choice < getNumChoices(state); choice++) {
			// Compute sum for this choice
			double d = mdpRewards.getTransitionReward(state, choice);
			d += mdpRewards.getStateReward(state);
			d += mvMultSingle(state, choice, vect);
			// Store strategy info if value matches
			if (PrismUtils.doublesAreCloseRel(val, d, PrismUtils.epsilonDouble)) {
				res.add(choice);
			}
		}

		return res;
	}

	/**
	 * Multiply the probability matrix induced by the MDP and {@code strat}
	 * to the right of {@code source}. Only those entries in {@code source}
	 * and only those columns in the probability matrix are considered, that
	 * are elements of {@code states}.
	 * <p>
	 * The result of this multiplication is added to the contents of {@code dest}.
	 *
	 * @param states States for which to multiply
	 * @param strat  (Memoryless) strategy to use
	 * @param source Vector to multiply matrix with
	 * @param dest   Vector to write result to.
	 */
	default void mvMultRight(int[] states, MDStrategy strat, double[] source, double[] dest)
	{
		for (int state : states) {
			for (Int2DoubleMap.Entry next : getTransitions(state, strat.getChoiceIndex(state))) {
				dest[next.getIntKey()] += next.getDoubleValue() * source[state];
			}
		}
	}

	@Override
	default String infoString()
	{
		String s = "";
		s += getNumStates() + " states (" + getNumInitialStates() + " initial)";
		s += ", " + getNumTransitions() + " transitions";
		s += ", " + getNumChoices() + " choices";
		s += ", dist max/avg = " + getMaxNumChoices() + "/" + PrismUtils.formatDouble2dp(((double) getNumChoices()) / getNumStates());
		return s;
	}

	@Override
	default String infoStringTable()
	{
		String s = "";
		s += "States:      " + getNumStates() + " (" + getNumInitialStates() + " initial)\n";
		s += "Transitions: " + getNumTransitions() + "\n";
		s += "Choices:     " + getNumChoices() + "\n";
		s += "Max/avg:     " + getMaxNumChoices() + "/" + PrismUtils.formatDouble2dp(((double) getNumChoices()) / getNumStates()) + "\n";
		return s;
	}

	@Override
	default Model constructInducedModel(MDStrategy strat)
	{
		return new DTMCFromMDPAndMDStrategy(this, strat);
	}

	@Override
	default void exportToPrismLanguage(String filename) throws PrismException
	{
		int i, j, numChoices;
		boolean first;
		Int2DoubleSortedMap sorted;
		Object action;
		int numStates = getNumStates();
		try (FileWriter out = new FileWriter(filename)) {
			// Output transitions to PRISM language file
			out.write(getModelType().keyword() + "\n");
			out.write("module M\nx : [0.." + (numStates - 1) + "];\n");
			sorted = new Int2DoubleAVLTreeMap();
			for (i = 0; i < numStates; i++) {
				numChoices = getNumChoices(i);
				for (j = 0; j < numChoices; j++) {
					// Extract transitions and sort by destination state index (to match PRISM-exported files)
					for (Int2DoubleMap.Entry e : getTransitions(i, j)) {
						sorted.put(e.getIntKey(), e.getDoubleValue());
					}
					// Print out (sorted) transitions
					action = getAction(i, j);
					out.write(action != null ? ("[" + action + "]") : "[]");
					out.write("x=" + i + "->");
					first = true;
					for (Int2DoubleMap.Entry e : sorted.int2DoubleEntrySet()) {
						if (first)
							first = false;
						else
							out.write("+");
						// Note use of PrismUtils.formatDouble to match PRISM-exported files
						out.write(PrismUtils.formatDouble(e.getDoubleValue()) + ":(x'=" + e.getIntKey() + ")");
					}
					out.write(";\n");
					sorted.clear();
				}
			}
			out.write("endmodule\n");
		} catch (IOException e) {
			throw new PrismException("Could not export " + getModelType() + " to file \"" + filename + "\"" + e);
		}
	}

	@Override
	default void exportToPrismExplicitTra(PrismLog out)
	{
		int i, j, numChoices;
		Object action;
		Int2DoubleSortedMap sorted;
		int numStates = getNumStates();
		// Output transitions to .tra file
		out.print(numStates + " " + getNumChoices() + " " + getNumTransitions() + "\n");
		sorted = new Int2DoubleAVLTreeMap();
		for (i = 0; i < numStates; i++) {
			numChoices = getNumChoices(i);
			for (j = 0; j < numChoices; j++) {
				// Extract transitions and sort by destination state index (to match PRISM-exported files)
				getTransitions(i, j).forEach(e -> sorted.put(e.getIntKey(), e.getDoubleValue()));

				// Print out (sorted) transitions
				for (Int2DoubleMap.Entry e : sorted.int2DoubleEntrySet()) {
					// Note use of PrismUtils.formatDouble to match PRISM-exported files
					out.print(i + " " + j + " " + e.getIntKey() + " " + PrismUtils.formatDouble(e.getDoubleValue()));
					action = getAction(i, j);
					out.print(action == null ? "\n" : (" " + action + "\n"));
				}
				sorted.clear();
			}
		}
	}

	@Override
	default void exportToDotFileWithStrat(PrismLog out, IntSet mark, MDStrategy strat)
	{
		int i, j, numChoices;
		String nij;
		Object action;
		String style;
		out.print("digraph " + getModelType() + " {\nsize=\"8,5\"\nnode [shape=box];\n");
		for (i = 0; i < getNumStates(); i++) {
			if (mark != null && mark.contains(i))
				out.print(i + " [style=filled  fillcolor=\"#cccccc\"]\n");
			numChoices = getNumChoices(i);
			for (j = 0; j < numChoices; j++) {
				style = (strat.getChoiceIndex(i) == j) ? ",color=\"#ff0000\",fontcolor=\"#ff0000\"" : "";
				action = getAction(i, j);
				nij = "n" + i + "_" + j;
				out.print(i + " -> " + nij + " [ arrowhead=none,label=\"" + j);
				if (action != null)
					out.print(":" + action);
				out.print("\"" + style + " ];\n");
				out.print(nij + " [ shape=point,height=0.1,label=\"\"" + style + " ];\n");
				for (Int2DoubleMap.Entry e : getTransitions(i, j)) {
					out.print(nij + " -> " + e.getIntKey() + " [ label=\"" + e.getDoubleValue() + "\" ];\n");
				}
			}
		}
		out.print("}\n");
	}

	enum Info
	{
		NO_TRANSIENT_MEC
	}

	enum Type
	{
		UNKNOWN, RECURRENT, UNICHAIN, COMMUNICATING, WEAKLY_COMMUNICATING, MULTICHAIN
	}

	@FunctionalInterface interface StateChoiceValueFunction
	{
		double value(int state, int choice);
	}
}
