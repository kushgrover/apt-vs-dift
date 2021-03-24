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

package pta;

import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;

import java.util.ArrayList;
import java.util.List;

public class LocZoneSet
{
	private Object2IntMap<LocZone> set;
	private int indexOfLastAdd;

	public LocZoneSet()
	{
		set = new Object2IntOpenHashMap<>();
		indexOfLastAdd = -1;
	}

	public boolean add(int loc, Zone z)
	{
		return add(new LocZone(loc, z));
	}

	public boolean add(LocZone lz)
	{
		if (set.containsKey(lz)) {
			indexOfLastAdd = set.getInt(lz);
			return false;
		} else {
			indexOfLastAdd = set.size();
			set.put(lz, set.size());
			return true;
		}
	}

	public boolean contains(LocZone lz)
	{
		return set.containsKey(lz);
	}

	public int getIndexOfLastAdd()
	{
		return indexOfLastAdd;
	}

	public boolean isEmpty()
	{
		return set.isEmpty();
	}

	public int size()
	{
		return set.size();
	}

	public List<LocZone> toList()
	{
		int i, n;

		n = set.size();
		List<LocZone> list = new ArrayList<>(n);
		for (i = 0; i < n; i++)
			list.add(null);
		for (Object2IntMap.Entry<LocZone> e : set.object2IntEntrySet()) {
			list.set(e.getIntValue(), e.getKey());
		}
		return list;
	}

	public String toString()
	{
		String s = "";
		s += set;
		return s;
	}
}
