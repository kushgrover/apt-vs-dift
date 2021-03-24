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

package heuristics.nextstate;

import heuristics.nextstate.HeuristicNextStateFactory.NextState;
import heuristics.update.StateUpdate;

import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import explicit.ModelExplorer;

import parser.State;
import prism.PrismException;


public class HeuristicNextStateDeterministic extends HeuristicNextState{

	private Map<State, Map<Integer,Map<Integer,Integer>>> state2LastTransBounded = new HashMap<State,Map<Integer,Map<Integer,Integer>>>();
	private Map<State, Map<Integer,Integer>> state2LastTrans = new HashMap<State,Map<Integer,Integer>>();
	
	
	public HeuristicNextStateDeterministic(ModelExplorer pme, StateUpdate stateUpdate) throws PrismException{
		super(pme, stateUpdate);
	}
	
	public State sample(State s, int action, int depth) throws PrismException{
		pme.queryState(s);
		int lastTrans = 0;
		if(state2LastTransBounded.get(s) != null) {
			if(state2LastTransBounded.get(s).get(depth) != null) {
				if(state2LastTransBounded.get(s).get(depth).get(action) != null) {
					lastTrans = state2LastTransBounded.get(s).get(depth).get(action);
					lastTrans++;
					int numTrans = pme.getNumTransitions(action);
					if(lastTrans >= numTrans) {
						lastTrans = 0;
					}
				}
			} else {
				state2LastTransBounded.get(s).put(depth,new HashMap<Integer,Integer>());
			}
		} else {
			state2LastTransBounded.put(s, new HashMap<Integer,Map<Integer,Integer>>());
			state2LastTransBounded.get(s).put(depth,new HashMap<Integer,Integer>());
		}
		state2LastTransBounded.get(s).get(depth).put(action, lastTrans);
		return pme.computeTransitionTarget(action, lastTrans);
	}
	
	public State sample(State s, int action) throws PrismException{
		pme.queryState(s);
		int lastTrans = 0;
		if(state2LastTrans.get(s) != null) {
			if(state2LastTrans.get(s).get(action) != null) {
				lastTrans = state2LastTrans.get(s).get(action);
				lastTrans++;
				int numTrans = pme.getNumTransitions(action);
				if(lastTrans >= numTrans) {
					lastTrans = 0;
				}
			}
		} else {
			state2LastTrans.put(s, new HashMap<Integer,Integer>());
		}
		state2LastTrans.get(s).put(action, lastTrans);
		return pme.computeTransitionTarget(action, lastTrans);
	}
	
	public NextState getType() {
		return NextState.DETERMINISTIC;
	}
				
}
