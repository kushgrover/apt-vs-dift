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

import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;

import parser.State;
import heuristics.search.StateValue;

public class StateValueContainer {
	
	protected Map<State,Map<Integer, StateValue>> qValues = new ConcurrentHashMap<State, Map<Integer,StateValue>>();
	protected Map<State,Integer> currentAction = new ConcurrentHashMap<State, Integer>();
	
	public void setQValue(State s, StateValue val, int depth) {
		if(qValues.get(s) == null) {
			Map<Integer, StateValue> depth2StateValue = new ConcurrentHashMap<Integer, StateValue>();
			depth2StateValue.put(depth, val);
			qValues.put(s, depth2StateValue);
		} else {
			qValues.get(s).put(depth, val);
		}
	}
	
	public StateValue getQValue(State s, int depth) {
		if(qValues.get(s) != null) {
			if(qValues.get(s).get(depth) != null) {
				return qValues.get(s).get(depth);
			} else {
				return null;
			}
		} else {
			return null;
		}
	}
	
	public void setCurrentAction(State s, int action) {
		currentAction.put(s, action);
	}
	
	public Integer getCurrentAction(State s) {
		return currentAction.get(s);
	}
	
	public int getSize() {
		return qValues.size();
	}
	
	public int getNumAllValues() {
		Iterator<Entry<State,Map<Integer, StateValue>>> it = qValues.entrySet().iterator();
		int allValues = 0;
		while(it.hasNext()) {
			Entry<State,Map<Integer, StateValue>> e = it.next();
			allValues += e.getValue().size();
		}
		return allValues;
	}
	
}
