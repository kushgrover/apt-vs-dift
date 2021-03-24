//==============================================================================
//
//	Copyright (c) 2013-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
//	* Mateusz Ujma <mateusz.ujma@cs.ox.ac.uk> (University of Oxford)
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

package heuristics.update;

import heuristics.search.StateValue;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import parser.State;

import java.util.HashMap;
import java.util.Map;

public final class DefaultDepthStateValueContainer implements StateValueContainer
{
	private final Map<State, Int2ObjectMap<StateValue>> qValues = new HashMap<>();
	private final Object2IntMap<State> currentAction = new Object2IntOpenHashMap<>();

	@Override public boolean setQValue(State state, StateValue value, int depth)
	{
		return qValues.computeIfAbsent(state, key -> new Int2ObjectOpenHashMap<>()).put(depth, value) != null;
	}

	@Override public StateValue getQValue(State state, int depth)
	{
		Int2ObjectMap<StateValue> valueMap = qValues.get(state);
		if (valueMap == null) {
			return null;
		}
		return valueMap.get(depth);
	}

	@Override public void setCurrentAction(State state, int action)
	{
		currentAction.put(state, action);
	}

	@Override public int getCurrentAction(State state)
	{
		return currentAction.getInt(state);
	}

	@Override public int getSize()
	{
		return qValues.size();
	}

	@Override public int getNumAllValues()
	{
		return qValues.values().stream().mapToInt(Map::size).sum();
	}
}
