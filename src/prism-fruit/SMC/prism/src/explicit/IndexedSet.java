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

import it.unimi.dsi.fastutil.objects.Object2IntAVLTreeMap;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import static com.google.common.base.Preconditions.checkArgument;

/**
 * Class storing an indexed set of objects of type T.
 * Typically used for storing state space during reachability.
 */
public class IndexedSet<T> implements StateStorage<T>
{
	protected Object2IntMap<T> set;
	protected int indexOfLastAdd;

	public IndexedSet()
	{
		this(false);
	}

	public IndexedSet(boolean sorted)
	{
		indexOfLastAdd = -1;
		set = sorted ? new Object2IntAVLTreeMap<>() : new Object2IntOpenHashMap<>();
		set.defaultReturnValue(-1);
	}

	@Override
	public void clear()
	{
		set.clear();
	}

	@Override
	public boolean add(T state)
	{
		int i = set.getInt(state);
		if (i > -1) {
			indexOfLastAdd = i;
			return false;
		} else {
			indexOfLastAdd = set.size();
			set.put(state, set.size());
			return true;
		}
	}

	@Override
	public boolean contains(T state)
	{
		return set.containsKey(state);
	}

	@Override
	public int getIndexOfLastAdd()
	{
		return indexOfLastAdd;
	}

	@Override
	public boolean isEmpty()
	{
		return set.isEmpty();
	}

	/**
	 * Get the number of objects stored in the set.
	 */
	@Override
	public int size()
	{
		return set.size();
	}

	/**
	 * Get access to the underlying set of map entries.
	 */
	@Override
	public Set<Object2IntMap.Entry<T>> getEntrySet()
	{
		return set.object2IntEntrySet();
	}

	/**
	 * Create an ArrayList of the states, ordered by index.
	 */
	@Override
	public List<T> toList()
	{
		List<T> list = new ArrayList<>(set.size());
		toList(list);
		return list;
	}

	/**
	 * Create a list of the states, ordered by index, storing in the passed in list.
	 *
	 * @param list An empty list in which to store the result.
	 */
	@Override
	public void toList(List<T> list)
	{
		checkArgument(list.isEmpty());

		int n = set.size();
		for (int i = 0; i < n; i++) {
			list.add(null);
		}
		for (Object2IntMap.Entry<T> e : set.object2IntEntrySet()) {
			list.set(e.getIntValue(), e.getKey());
		}
	}

	/**
	 * Create an ArrayList of the states, ordered by permuted index.
	 * Index in new list is permut[old_index].
	 *
	 * @param permut Permutation to apply
	 */
	@Override
	public List<T> toPermutedList(int permut[])
	{
		List<T> list = new ArrayList<>(set.size());
		toPermutedList(permut, list);
		return list;
	}

	/**
	 * Create an ArrayList of the states, ordered by permuted index, storing in the passed in list.
	 * Index in new list is permut[old_index].
	 *
	 * @param permut Permutation to apply
	 * @param list   An empty ArrayList in which to store the result.
	 */
	@Override
	public void toPermutedList(int permut[], List<T> list)
	{
		int i, n;

		n = set.size();
		for (i = 0; i < n; i++) {
			list.add(null);
		}
		for (Object2IntMap.Entry<T> e : set.object2IntEntrySet()) {
			list.set(permut[e.getIntValue()], e.getKey());
		}
	}

	/**
	 * Build sort permutation. Assuming this was built as a sorted set,
	 * this returns a permutation (integer array) mapping current indices
	 * to new indices under the sorting order.
	 */
	@Override
	public int[] buildSortingPermutation()
	{
		int i, n;
		int perm[];

		n = set.size();
		perm = new int[n];
		i = 0;
		for (Object2IntMap.Entry<T> e : set.object2IntEntrySet()) {
			perm[e.getIntValue()] = i;
			i++;
		}

		return perm;
	}

	@Override
	public String toString()
	{
		return set.toString();
	}

	@Override
	public int get(T t)
	{
		return set.getInt(t);
	}
}
