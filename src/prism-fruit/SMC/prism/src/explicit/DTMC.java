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

import com.google.common.collect.Iterators;
import common.FastUtils;
import de.tum.in.naturals.set.NatBitSet;
import it.unimi.dsi.fastutil.ints.AbstractInt2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleAVLTreeMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleSortedMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectAVLTreeMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectSortedMap;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntSet;
import prism.ModelType;
import prism.Pair;
import prism.PrismException;
import prism.PrismLog;
import prism.PrismUtils;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.function.IntConsumer;
import java.util.function.IntToDoubleFunction;

import static com.google.common.base.Preconditions.checkArgument;

/**
 * Interface for classes that provide (read) access to an explicit-state DTMC.
 */
public interface DTMC extends Model
{
	/**
	 * Get the number of transitions from state s.
	 */
	default int getNumTransitions(int s) {
		return Iterators.size(getTransitionsIterator(s));
	}

	@Override
	default void exportToPrismExplicitTra(PrismLog out)
	{
		int numStates = getNumStates();
		Int2ObjectSortedMap<Pair<Double, Object>> sorted = new Int2ObjectAVLTreeMap<>();
		// Output transitions to .tra file
		out.print(numStates + " " + getNumTransitions() + "\n");
		for (int i = 0; i < numStates; i++) {
			// Extract transitions and sort by destination state index (to match PRISM-exported files)
			getTransitionsAndActionsIterator(i).forEachRemaining(e -> sorted.put(e.getIntKey(), e.getValue()));
			// Print out (sorted) transitions
			for (Int2ObjectMap.Entry<Pair<Double, Object>> e : sorted.int2ObjectEntrySet()) {
				// Note use of PrismUtils.formatDouble to match PRISM-exported files
				out.print(i + " " + e.getIntKey() + " " + PrismUtils.formatDouble(e.getValue().first));
				Object action = e.getValue().second;
				if (action != null && !"".equals(action)) {
					out.print(" " + action);
				}
				out.print("\n");
			}
			sorted.clear();
		}
	}

	@Override
	default void exportToPrismLanguage(String filename) throws PrismException
	{
		int numStates = getNumStates();
		try (FileWriter out = new FileWriter(filename)) {
			// Output transitions to PRISM language file
			out.write(getModelType().keyword() + "\n");
			out.write("module M\nx : [0.." + (numStates - 1) + "];\n");
			Int2DoubleSortedMap sorted = new Int2DoubleAVLTreeMap();
			for (int i = 0; i < numStates; i++) {
				// Extract transitions and sort by destination state index (to match PRISM-exported files)
				Iterator<Int2DoubleMap.Entry> iter = getTransitionsIterator(i);
				while (iter.hasNext()) {
					Int2DoubleMap.Entry e = iter.next();
					sorted.put(e.getIntKey(), e.getDoubleValue());
				}
				// Print out (sorted) transitions
				out.write("[]x=" + i + "->");
				boolean first = true;
				for (Int2DoubleMap.Entry e : sorted.int2DoubleEntrySet()) {
					if (first) {
						first = false;
					} else {
						out.write("+");
					}
					// Note use of PrismUtils.formatDouble to match PRISM-exported files
					out.write(PrismUtils.formatDouble(e.getDoubleValue()) + ":(x'=" + e.getIntKey() + ")");
				}
				out.write(";\n");
				sorted.clear();
			}
			out.write("endmodule\n");
			out.close();
		} catch (IOException e) {
			throw new PrismException("Could not export " + getModelType() + " to file \"" + filename + "\"" + e);
		}
	}

	/**
	 * Get an iterator over the transitions from state s.
	 */
	Iterator<Int2DoubleMap.Entry> getTransitionsIterator(int s);

	/**
	 * Get an iterator over the transitions from state s, with their attached actions if present.
	 */
	default Iterator<Int2ObjectMap.Entry<Pair<Double, Object>>> getTransitionsAndActionsIterator(int s)
	{
		return new AddDefaultActionToTransitionsIterator(getTransitionsIterator(s), null);
	}

	/**
	 * Perform a single step of precomputation algorithm Prob0, i.e., for states i in {@code subset},
	 * set bit i of {@code result} iff there is a transition to a state in {@code u}.
	 *
	 * @param subset Only compute for these states
	 * @param u      Set of states {@code u}
	 * @param result Store results here
	 */
	default void prob0step(IntSet subset, IntSet u, NatBitSet result)
	{
		subset.forEach((IntConsumer) i -> {
			IntIterator successorsIterator = getSuccessorsIterator(i);
			while (successorsIterator.hasNext()) {
				if (u.contains(successorsIterator.nextInt())) {
					result.set(i);
					return;
				}
			}
			result.clear(i);
		});
	}

	/**
	 * Perform a single step of precomputation algorithm Prob1, i.e., for states i in {@code subset},
	 * set bit i of {@code result} iff there is a transition to a state in {@code v} and all transitions go to states in {@code u}.
	 *
	 * @param subset Only compute for these states
	 * @param u      Set of states {@code u}
	 * @param v      Set of states {@code v}
	 * @param result Store results here
	 */
	default void prob1step(IntSet subset, IntSet u, IntSet v, NatBitSet result)
	{
		subset.forEach((IntConsumer) i -> {
			IntIterator successorsIterator = getSuccessorsIterator(i);
			boolean hasV = false;
			while (successorsIterator.hasNext()) {
				int successor = successorsIterator.nextInt();
				if (!(u.contains(successor))) {
					result.clear(i);
					return;
				}
				hasV |= v.contains(successor);
			}
			result.set(i, hasV);
		});
	}

	/**
	 * Do a matrix-vector multiplication for
	 * the DTMC's transition probability matrix P and the vector {@code vect} passed in.
	 * i.e. for all s: result[s] = sum_j P(s,j)*vect[j]
	 *
	 * @param vect   Vector to multiply by
	 * @param result Vector to store result in
	 * @param subset Only do multiplication for these rows
	 */
	default void mvMult(double vect[], double result[], IntSet subset)
	{
		FastUtils.forEach(subset, getNumStates(), s -> result[s] = mvMultSingle(s, vect));
	}

	/**
	 * Do a single row of matrix-vector multiplication for
	 * the DTMC's transition probability matrix P and the vector {@code vect} passed in.
	 * i.e. return sum_j P(s,j)*vect[j]
	 *
	 * @param s    Row index
	 * @param vect Vector to multiply by
	 */
	default double mvMultSingle(int s, IntToDoubleFunction vect)
	{
		Iterator<Int2DoubleMap.Entry> transitionsIterator = getTransitionsIterator(s);
		double sum = 0.0d;
		while (transitionsIterator.hasNext()) {
			Int2DoubleMap.Entry transition = transitionsIterator.next();
			sum += transition.getDoubleValue() * vect.applyAsDouble(transition.getIntKey());
		}
		return sum;
	}

	@Override
	default ModelType getModelType()
	{
		return ModelType.DTMC;
	}

	/**
	 * Do a single row of matrix-vector multiplication for
	 * the DTMC's transition probability matrix P and the vector {@code vect} passed in.
	 * i.e. return sum_j P(s,j)*vect[j]
	 *
	 * @param s    Row index
	 * @param vect Vector to multiply by
	 */
	default double mvMultSingle(int s, double[] vect)
	{
		return mvMultSingle(s, state -> vect[state]);
	}

	/**
	 * Do a Gauss-Seidel-style matrix-vector multiplication for
	 * the DTMC's transition probability matrix P and the vector {@code vect} passed in,
	 * storing new values directly in {@code vect} as computed.
	 * i.e. for all s: vect[s] = (sum_{j!=s} P(s,j)*vect[j]) / (1-P(s,s))
	 * The maximum (absolute/relative) difference between old/new
	 * elements of {@code vect} is also returned.
	 *
	 * @param vect     Vector to multiply by (and store the result in)
	 * @param subset   Only do multiplication for these rows (ignored if null)
	 * @param absolute If true, compute absolute, rather than relative, difference
	 * @return The maximum difference between old/new elements of {@code vect}
	 */
	default double mvMultGS(double vect[], IntSet subset, boolean absolute)
	{
		double d, diff, maxDiff = 0.0;
		IntIterator iterator = subset.iterator();
		while (iterator.hasNext()) {
			int s = iterator.nextInt();
			d = mvMultJacSingle(s, vect);
			diff = absolute ? (Math.abs(d - vect[s])) : (Math.abs(d - vect[s]) / d);
			maxDiff = diff > maxDiff ? diff : maxDiff;
			vect[s] = d;
		}
		// Use this code instead for backwards Gauss-Seidel
		/*for (s = numStates - 1; s >= 0; s--) {
			if (subset.get(s)) {
				d = mvMultJacSingle(s, vect);
				diff = absolute ? (Math.abs(d - vect[s])) : (Math.abs(d - vect[s]) / d);
				maxDiff = diff > maxDiff ? diff : maxDiff;
				vect[s] = d;
			}
		}*/
		return maxDiff;
	}

	/**
	 * Do a single row of Jacobi-style matrix-vector multiplication for
	 * the DTMC's transition probability matrix P and the vector {@code vect} passed in.
	 * i.e. return (sum_{j!=s} P(s,j)*vect[j]) / (1-P(s,s))
	 *
	 * @param s    Row index
	 * @param vect Vector to multiply by
	 */
	default double mvMultJacSingle(int s, double vect[])
	{
		double diag = 1.0d;
		double d = 0.0d;
		Iterator<Int2DoubleMap.Entry> transitionsIterator = getTransitionsIterator(s);

		while (transitionsIterator.hasNext()) {
			Int2DoubleMap.Entry transition = transitionsIterator.next();
			int successor = transition.getIntKey();
			double prob = transition.getDoubleValue();
			if (successor != s) {
				d += prob * vect[successor];
			} else {
				diag -= prob;
			}
		}
		if (diag > 0) {
			d /= diag;
		}
		return d;
	}

	/**
	 * Performs {@link #mvMultSingle(int, IntToDoubleFunction)} on DTMC after performing aperiodicity transformation
	 * as per Section 8.5.4 of Puterman
	 */
	default double mvMultSingleAperiodic(int state, IntToDoubleFunction values, double tau)
	{
		checkArgument(0 < tau && tau <= 1d, "Aperiodicity parameter has to be in (0, 1]");
		double successorValues = mvMultSingle(state, values);
		if (tau == 1.0d) {
			return successorValues;
		}
		return (1 - tau) * values.applyAsDouble(state) + tau * successorValues;
	}

	default double mvMultSingleAperiodic(int state, double values[], double tau)
	{
		return mvMultSingleAperiodic(state, s -> values[s], tau);
	}

	/**
	 * Do a vector-matrix multiplication for
	 * the DTMC's transition probability matrix P and the vector {@code vect} passed in.
	 * i.e. for all s: result[s] = sum_i P(i,s)*vect[i]
	 *
	 * @param vect   Vector to multiply by
	 * @param result Vector to store result in
	 */
	default void vmMult(double vect[], double result[])
	{
		// Initialise result to 0
		Arrays.fill(result, 0d);
		// Go through matrix elements (by row)
		for (int i = 0; i < getNumStates(); i++) {
			Iterator<Int2DoubleMap.Entry> transitionsIterator = getTransitionsIterator(i);
			while (transitionsIterator.hasNext()) {
				Int2DoubleMap.Entry transition = transitionsIterator.next();
				result[transition.getIntKey()] += transition.getDoubleValue() * vect[i];
			}
		}
	}

	default Type getMCType()
	{
		return Type.UNKNOWN;
	}

	enum Type
	{
		STRONGLY_UNICHAIN, UNICHAIN, UNKNOWN
	}

	class AddDefaultActionToTransitionsIterator implements Iterator<Int2ObjectMap.Entry<Pair<Double, Object>>>
	{
		private final Iterator<Int2DoubleMap.Entry> transIter;
		private final Object defaultAction;
		private Int2DoubleMap.Entry next;

		public AddDefaultActionToTransitionsIterator(Iterator<Int2DoubleMap.Entry> transIter, Object defaultAction)
		{
			this.transIter = transIter;
			this.defaultAction = defaultAction;
		}

		@Override
		public Int2ObjectMap.Entry<Pair<Double, Object>> next()
		{
			next = transIter.next();
			int state = next.getIntKey();
			double probability = next.getDoubleValue();
			return new AbstractInt2ObjectMap.BasicEntry<>(state, new Pair<>(probability, defaultAction));
		}

		@Override
		public boolean hasNext()
		{
			return transIter.hasNext();
		}
	}
}
