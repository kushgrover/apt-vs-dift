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

import com.google.common.collect.Ordering;
import common.IntDoubleConsumer;
import it.unimi.dsi.fastutil.ints.Int2DoubleFunction;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntCollection;
import prism.PrismUtils;

import java.util.Arrays;
import java.util.Iterator;
import java.util.function.IntConsumer;
import java.util.function.IntToDoubleFunction;
import java.util.function.IntUnaryOperator;

import static com.google.common.base.Preconditions.checkArgument;

/**
 * Explicit representation of a probability distribution.
 * Basically, a mapping from (integer-valued) indices to (non-zero, double-valued) probabilities.
 */
public final class Distribution implements Iterable<Int2DoubleMap.Entry>
{
	private static final int[] INTS = new int[] {};
	private static final double[] DOUBLES = new double[] {};

	private int[] keys;
	private double[] values;

	public static Distribution create() {
		return new Distribution();
	}

	public static Distribution create(int key, double value) {
		return new Distribution(key, value);
	}

	public static Distribution copy(Distribution distribution) {
		return new Distribution(distribution);
	}

	public static Distribution permute(Distribution distribution, int[] permutation) {
		DistributionBuilder build = new DistributionBuilder(distribution.size());
		for (int i = 0, keysLength = distribution.keys.length; i < keysLength; i++) {
			build.add(permutation[distribution.keys[i]], distribution.values[i]);
		}
		return build.build();
	}

	/**
	 * Create an empty distribution.
	 */
	public Distribution()
	{
		keys = INTS;
		values = DOUBLES;
	}

	public Distribution(int key, double value)
	{
		keys = new int[] { key };
		values = new double[] { value };
	}

	/**
	 * Copy constructor.
	 */
	public Distribution(Distribution distr)
	{
		keys = distr.keys.clone();
		values = distr.values.clone();
	}

	/**
	 * Construct a distribution from an existing one and an index permutation,
	 * i.e. in which index i becomes index permut[i].
	 * Note: have to build the new distributions from scratch anyway to do this,
	 * so may as well provide this functionality as a constructor.
	 */
	public Distribution(Distribution distr, int permut[])
	{
		keys = INTS;
		values = DOUBLES;
		distr.forEach((key, value) -> add(permut[key], value));
	}

	private Distribution(int[] keys, double[] values)
	{
		assert keys.length == values.length && Ordering.natural().isStrictlyOrdered(new IntArrayList(keys));
		this.keys = keys;
		this.values = values;
	}

	/**
	 * Clear all entries of the distribution.
	 */
	public void clear()
	{
		keys = INTS;
		values = DOUBLES;
	}

	/**
	 * Add 'prob' to the probability for index 'j'.
	 * Return boolean indicating whether or not there was already
	 * non-zero probability for this index (i.e. false denotes new transition).
	 *
	 * @param j    - the state j
	 * @param prob - the probability mass added to state j.
	 * @return
	 */
	public boolean add(int j, double prob)
	{
		checkArgument(j >= 0, "State index is negative: %s", j);
		checkArgument(prob >= 0.d, "Probability is negative: %s", prob);
		if (prob == 0.d) {
			// No new transition
			return true;
		}

		int i = Arrays.binarySearch(keys, j);

		if (i >= 0) {
			// Have this key already
			values[i] += prob;
			return true;
		}
		// Need to insert it

		// Compute insertion point
		i = -(i + 1);

		int[] localKeys = new int[keys.length + 1];
		double[] localValues = new double[keys.length + 1];

		System.arraycopy(keys, 0, localKeys, 0, i);
		System.arraycopy(values, 0, localValues, 0, i);
		localKeys[i] = j;
		localValues[i] = prob;
		System.arraycopy(keys, i, localKeys, i + 1, keys.length - i);
		System.arraycopy(values, i, localValues, i + 1, keys.length - i);

		keys = localKeys;
		values = localValues;

		return false;
	}

	/**
	 * Set the probability for index 'j' to 'prob'.
	 */
	public void set(int j, double prob)
	{
		checkArgument(prob > 0, "Probability is non-positive: %s", prob);

		int i = Arrays.binarySearch(keys, j);

		if (i >= 0) {
			values[i] = prob;
		} else {
			add(j, prob);
		}
	}

	/**
	 * Get the probability for index j.
	 */
	public double get(int j)
	{
		int position = Arrays.binarySearch(keys, j);

		if (position < 0) {
			return 0.;
		}

		return values[position];
	}

	/**
	 * Returns true if index j is in the support of the distribution.
	 */
	public boolean contains(int j)
	{
		return Arrays.binarySearch(keys, j) >= 0;
	}

	/**
	 * Returns true if all indices in the support of the distribution are in the set.
	 */
	public boolean isSubsetOf(IntCollection set)
	{
		for (int key : keys) {
			if (!set.contains(key)) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Returns true if at least one index in the support of the distribution is in the set.
	 */
	public boolean containsOneOf(IntCollection set)
	{
		for (int key : keys) {
			if (set.contains(key)) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Get the support of the distribution.
	 */
	public int[] getSupport()
	{
		return keys;
	}

	/**
	 * Get an iterator over the entries of the map defining the distribution.
	 */
	@Override public Iterator<Int2DoubleMap.Entry> iterator()
	{
		return new Iterator<Int2DoubleMap.Entry>()
		{
			private final DistributionEntry entry = new DistributionEntry();

			@Override
			public boolean hasNext()
			{
				return (entry.index + 1) < keys.length;
			}

			@Override
			public Int2DoubleMap.Entry next()
			{
				entry.index++;
				return entry;
			}
		};
	}

	/**
	 * Returns true if the distribution is empty.
	 */
	public boolean isEmpty()
	{
		return keys.length == 0;
	}

	/**
	 * Get the size of the support of the distribution.
	 */
	public int size()
	{
		return keys.length;
	}

	/**
	 * Get the sum of the probabilities in the distribution.
	 */
	public double sum()
	{
		double d = 0.0d;
		for (double value : values) {
			d += value;
		}
		return d;
	}

	public double sumWeighted(double[] weights)
	{
		double d = 0.0d;
		for (int i = 0, keysLength = keys.length; i < keysLength; i++) {
			d += weights[keys[i]] * values[i];
		}
		return d;
	}

	public double sumWeighted(Int2DoubleFunction weights)
	{
		double d = 0.0d;
		for (int i = 0, keysLength = keys.length; i < keysLength; i++) {
			d += weights.get(keys[i]) * values[i];
		}
		return d;
	}

	public double sumWeighted(IntToDoubleFunction weights)
	{
		double d = 0.0d;
		for (int i = 0, keysLength = keys.length; i < keysLength; i++) {
			d += weights.applyAsDouble(keys[i]) * values[i];
		}
		return d;
	}

	/**
	 * Get the sum of all the probabilities in the distribution except for index j.
	 */
	public double sumAllBut(int j)
	{
		double d = 0.0d;
		for (int i = 0, keysLength = keys.length; i < keysLength; i++) {
			if (keys[i] != j) {
				d += values[i];
			}
		}
		return d;
	}

	/**
	 * Create a new distribution, based on a mapping from the indices
	 * used in this distribution to a different set of indices.
	 */
	public Distribution map(int map[])
	{
		DistributionBuilder builder = new DistributionBuilder(size());
		for (int i = 0, keysLength = keys.length; i < keysLength; i++) {
			builder.add(map[keys[i]], values[i]);
		}
		return builder.build();
	}

	/**
	 * Create a new distribution, based on a mapping from the indices
	 * used in this distribution to a different set of indices.
	 */
	public Distribution map(IntUnaryOperator map)
	{
		DistributionBuilder builder = new DistributionBuilder(size());
		for (int i = 0, keysLength = keys.length; i < keysLength; i++) {
			builder.add(map.applyAsInt(keys[i]), values[i]);
		}
		return builder.build();
	}

	@Override
	public boolean equals(Object o)
	{
		if (this == o)
			return true;
		if (o == null || getClass() != o.getClass())
			return false;
		Distribution that = (Distribution) o;
		return Arrays.equals(keys, that.keys) &&
				PrismUtils.doublesAreClose(values, that.values, 1e-12, false);
	}

	@Override
	public int hashCode()
	{
		return Arrays.hashCode(keys);
	}

	@Override
	public String toString()
	{
		StringBuilder builder = new StringBuilder();
		builder.append('{');
		for (Iterator<Int2DoubleMap.Entry> iterator = this.iterator(); iterator.hasNext(); ) {
			Int2DoubleMap.Entry entry = iterator.next();
			builder.append(entry.getIntKey()).append('=').append(entry.getDoubleValue());
			if (iterator.hasNext()) {
				builder.append(", ");
			}
		}
		builder.append('}');

		return builder.toString();
	}

	public void forEach(IntDoubleConsumer consumer)
	{
		for (int i = 0, keysLength = keys.length; i < keysLength; i++) {
			consumer.accept(keys[i], values[i]);
		}
	}

	public void forEach(IntConsumer consumer)
	{
		for (int key : keys) {
			consumer.accept(key);
		}
	}

	public void scale(double scale)
	{
		checkArgument(scale > 0d);
		for (int i = 0; i < values.length; i++) {
			values[i] *= scale;
		}
	}

	private static final class DistributionBuilder
	{
		private static final int INITIAL_SIZE = 16;
		private int[] keys;
		private double[] values;

		DistributionBuilder()
		{
			this(INITIAL_SIZE);
		}

		DistributionBuilder(int size)
		{
			keys = new int[size];
			values = new double[size];
			Arrays.fill(keys, Integer.MAX_VALUE);
		}

		private boolean add(int j, double prob)
		{
			checkArgument(j >= 0, "State index is negative: %s", j);
			checkArgument(prob > 0, "Probability is non-positive: %s", prob);

			int i = Arrays.binarySearch(keys, j);

			if (i >= 0) {
				// Have this key already
				values[i] += prob;
				return true;
			}

			i = -(i + 1);

			if (keys[keys.length - 1] != Integer.MAX_VALUE) {
				// No space left, resize
				int oldSize = keys.length;
				int newSize = (int) (oldSize * 1.5) + 1;
				keys = Arrays.copyOf(keys, newSize);
				values = Arrays.copyOf(values, newSize);
				Arrays.fill(keys, oldSize, newSize, Integer.MAX_VALUE);
			} else if (keys[i] == Integer.MAX_VALUE) {
				// Insertion point is free, overwrite it
				assert i == 0 || keys[i - 1] < Integer.MAX_VALUE;
				keys[i] = j;
				values[i] = prob;
				return false;
			}
			System.arraycopy(keys, i, keys, i + 1, keys.length - i - 1);
			System.arraycopy(values, i, values, i + 1, keys.length - i - 1);
			keys[i] = j;
			values[i] = prob;
			return false;
		}

		private Distribution build()
		{
			if (keys.length == 0 || keys[0] == Integer.MAX_VALUE) {
				// No entries at all
				return new Distribution();
			}
			if (keys[keys.length - 1] != Integer.MAX_VALUE) {
				// Our array is nicely filled
				return new Distribution(keys, values);
			}

			// Compact the array, removing superfluous indices
			assert keys.length >= 2;
			int i = keys.length - 2;
			while (keys[i] == Integer.MAX_VALUE) {
				i--;
			}
			int[] keys = Arrays.copyOf(this.keys, i + 1);
			double[] values = Arrays.copyOf(this.values, i + 1);
			return new Distribution(keys, values);
		}

	}

	private class DistributionEntry implements Int2DoubleMap.Entry
	{
		private int index = -1;

		@Override
		public int getIntKey()
		{
			return keys[index];
		}

		@Deprecated
		@Override public Integer getKey()
		{
			return getIntKey();
		}

		@Override
		public double setValue(double prob)
		{
			checkArgument(prob > 0., "Probability is out of range: " + prob);
			double oldValue = values[index];
			values[index] = prob;
			return oldValue;
		}

		@Deprecated
		@Override public Double getValue()
		{
			return getDoubleValue();
		}

		@Deprecated
		@Override public Double setValue(Double aDouble)
		{
			return setValue(aDouble.doubleValue());
		}

		@Override
		public double getDoubleValue()
		{
			return values[index];
		}
	}
}
