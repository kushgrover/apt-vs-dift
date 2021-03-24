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
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import parser.State;

import java.util.HashMap;
import java.util.Map;

public final class DefaultStateValueContainer implements StateValueContainer
{
	private final Map<State, StateValue> qValues = new HashMap<>();
	private final Object2IntMap<State> currentAction = new Object2IntOpenHashMap<>();

	@Override public boolean setQValue(State state, StateValue value)
	{
		return qValues.put(state, value) != null;
	}

	@Override public StateValue getQValue(State state)
	{
		return qValues.get(state);
	}

	@Override public StateValue getQValue(State state, int depth)
	{
		if (depth == -1) {
			return getQValue(state);
		}
		throw new UnsupportedOperationException(String.format("Depth unsupported (%d requested)", depth));
	}

	@Override public void setCurrentAction(State state, int action)
	{
		currentAction.put(state, action);
	}

	@Override public int getCurrentAction(State state)
	{
		return currentAction.getInt(state);
	}

	@Override public boolean setQValue(State state, StateValue value, int depth)
	{
		if (depth == -1) {
			return setQValue(state, value);
		}
		throw new UnsupportedOperationException(String.format("Depth unsupported (%d requested)", depth));
	}

	@Override public int getSize()
	{
		return qValues.size();
	}

	@Override public int getNumAllValues()
	{
		return qValues.size();
	}
}
