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

import java.util.*;

import explicit.STPGExplicit;
import heuristics.search.HeuristicRTDP_Adj;
import heuristics.search.StateValue;

import explicit.ModelExplorer;
import parser.State;
import prism.Pair;
import prism.PrismException;

public class StateUpdateMDP extends StateUpdate {

	public StateUpdateMDP(ModelExplorer me, StateValueContainer container, int bound, boolean min, double epsilon) throws PrismException {
		super(me, container, bound, min, epsilon);
	}

	private int randConstant = 1000000000;
	private int updateCount = 0;

	public int update(State s, int depth) throws PrismException {
		pme.queryState(s);
		int choices = pme.getNumChoices();
		StateValue lowerBound_q_s[] = new StateValue[choices];
		StateValue upperBound_q_s[] = new StateValue[choices];
		double lowerBound;
		double upperBound;
		if (min) {
			lowerBound = Double.MAX_VALUE;
			upperBound = Double.MAX_VALUE;
		} else {
			lowerBound = -1;
			upperBound = -1;
		}

		int bestAction = 0;
		for (int i = 0; i < choices; i++) {
			lowerBound_q_s[i] = getLowerBoundActionValue(s, i, depth + 1);
			upperBound_q_s[i] = getUpperBoundActionValue(s, i, depth + 1);
			if (min) {
				if (lowerBound_q_s[i].getLowerBound() < lowerBound) {
					bestAction = i;
					lowerBound = lowerBound_q_s[i].getLowerBound();
					upperBound = upperBound_q_s[i].getUpperBound();
				}
			} else {
				if (upperBound_q_s[i].getUpperBound() > upperBound) {
					bestAction = i;
					lowerBound = lowerBound_q_s[i].getLowerBound();
					upperBound = upperBound_q_s[i].getUpperBound();
				}
			}
		}

		setQValue(s, new StateValue(lowerBound, upperBound), depth);
		return bestAction;
	}

	public int update_Mat(State s) throws PrismException {
		pme.queryState(s);
		int choices = pme.getNumChoices();
		StateValue lowerBound_q_s[] = new StateValue[choices];
		StateValue upperBound_q_s[] = new StateValue[choices];
		double lowerBound;
		double upperBound;
		if (min) {
			lowerBound = Double.MAX_VALUE;
			upperBound = Double.MAX_VALUE;
		} else {
			lowerBound = -1;
			upperBound = -1;
		}
		int bestAction = 0;
		for (int i = 0; i < choices; i++) {
			lowerBound_q_s[i] = getLowerBoundActionValue(s, i);
			upperBound_q_s[i] = getUpperBoundActionValue(s, i);
			if (min) {
				if (lowerBound_q_s[i].getLowerBound() < lowerBound) {
					bestAction = i;
					lowerBound = lowerBound_q_s[i].getLowerBound();
					upperBound = upperBound_q_s[i].getUpperBound();
				}
			} else {
				if (upperBound_q_s[i].getUpperBound() > upperBound) {
					bestAction = i;
					lowerBound = lowerBound_q_s[i].getLowerBound();
					upperBound = upperBound_q_s[i].getUpperBound();
				}
			}
		}

		List<Integer> same = new ArrayList<Integer>();
		for (int i = 0; i < choices; i++) {
			lowerBound_q_s[i] = getLowerBoundActionValue(s, i);
			upperBound_q_s[i] = getUpperBoundActionValue(s, i);
			if (min) {
				if (same(lowerBound_q_s[i].getLowerBound(), lowerBound)) {
					same.add(i);
				}
			} else {
				if (same(upperBound_q_s[i].getUpperBound(), upperBound)) {
					same.add(i);
				}
			}
		}

		Random r = new Random();
		if (same.size() > 0) {
			bestAction = same.get(r.nextInt(same.size()));
		}

//		StateValue sv = getQValue(s);
//		if (sv!=null && sv.getLowerBound()>getLowerBoundActionValue(s, bestAction).getLowerBound()){
//			int breakpoint = 1;
//		}
		setQValue(s,
				new StateValue(getLowerBoundActionValue(s, bestAction).getLowerBound(),
						getUpperBoundActionValue(s, bestAction).getUpperBound()));

		updateCount++;
		return bestAction;
	}

	public StateUpdate clone() {
		try {
			return new StateUpdateMDP(null, null, bound, min, epsilon);
		} catch (PrismException e) {
			throw new RuntimeException(e);
		}
	}

	//A change in this needs to be changed in StateUpdateSMG as well
	public int update(State s, boolean isMaxState) throws PrismException {
		Pair<ArrayList<Integer>,StateValue> bestActionsAndValue = getBestActions(s,isMaxState);
		ArrayList<Integer> bestActions = bestActionsAndValue.first;
		//set QValue of the state according to what we just found out about best actions
		setQValue(s, bestActionsAndValue.second);
		//pick a random action from the set of best actions
		return bestActions.get(new Random().nextInt(bestActions.size()));
	}

	public Pair<ArrayList<Integer>,StateValue> getBestActions(State s, boolean isMaxState) throws PrismException {
		pme.queryState(s);
		int choices = pme.getNumChoices();
		ArrayList<Integer> bestActions = new ArrayList<Integer>();
		double currentLowerBound, currentUpperBound;
		double lowerBoundValue;
		double upperBoundValue;
		if (min) {
			lowerBoundValue = 1;
			upperBoundValue = 1;
		} else {
			lowerBoundValue = 0;
			upperBoundValue = 0;
		}
		for(int i=0;i<choices;i++) {
			pme.queryState(s);
			currentLowerBound = getLowerBoundActionValue(s, i).getLowerBound();
			currentUpperBound = getUpperBoundActionValue(s, i).getUpperBound();
			if (min) {
				if (currentUpperBound<upperBoundValue){
					upperBoundValue = currentUpperBound; //we always remember the minimal upper bound possible
				}
				if(lowerBoundValue - currentLowerBound > 0) { //the current action has a value that is smaller
					bestActions.clear();
					bestActions.add(i);
					lowerBoundValue = currentLowerBound;
				}
				else if (lowerBoundValue - currentLowerBound==0){//the current action is as good as best we have
					//another best Action, add it to the set
					bestActions.add(i);
				}
			} else {
				if (currentLowerBound>lowerBoundValue){
					lowerBoundValue = currentLowerBound; //we always remember the maximal lower bound possible
				}
				if(currentUpperBound - upperBoundValue>0) { //The current action is better by more than epsilon
					//new best Action, remove all others we had until now
					bestActions.clear();
					bestActions.add(i);
					upperBoundValue = currentUpperBound;
				}
				else if (currentUpperBound - upperBoundValue==0){ //The current action is same as best we have
					//another best Action, add it to the set
					bestActions.add(i);
				}
			}
		}
		assert(!bestActions.isEmpty());

		return new Pair<>(bestActions,new StateValue(lowerBoundValue, upperBoundValue));
	}

	//Relic. Maybe there was a reason, so I keep it.
	public ArrayList<Integer> getEpsBestActions(State s) throws PrismException {
		pme.queryState(s);
		int choices = pme.getNumChoices();
		ArrayList<Integer> bestActions = new ArrayList<Integer>();
		double currentLowerBound, currentUpperBound;
		double lowerBoundValue;
		double upperBoundValue;
		if(!min) {
			lowerBoundValue = 0;
			upperBoundValue = 0;
		} else {
			lowerBoundValue = 1;
			upperBoundValue = 1;
		}
		for(int i=0;i<choices;i++) {
			pme.queryState(s);
			currentLowerBound = getLowerBoundActionValue(s, i).getLowerBound();
			currentUpperBound = getUpperBoundActionValue(s, i).getUpperBound();
			if(!min) {
				if (currentLowerBound>lowerBoundValue){
					lowerBoundValue = currentLowerBound; //we always remember the maximal lower bound possible
				}
				if(currentUpperBound - upperBoundValue>epsilon) { //The current action is better by more than epsilon
					//new best Action, remove all others we had until now
					bestActions.clear();
					bestActions.add(i);
					upperBoundValue = currentUpperBound;
				}
				else if (Math.abs(currentUpperBound - upperBoundValue)<epsilon){ //The current action is epsilon close to the one that was best before
					//another best Action, add it to the set
					bestActions.add(i);
					//We want to stay epsilon close to the best action, so we need to remember the highest upperBound in bestActions, and we need to check that nothing forbidden is in bestActions
					if(currentUpperBound > upperBoundValue){
						upperBoundValue = currentUpperBound;
						lowerBoundValue = currentLowerBound; //we remember the lower bound of the action with the maximal U. This is definitely correct.
						//The new action is greater than the old; therefor the epsilon area of the best action shifts, and we have to check, whether we have to remove something from bestActions
						for (int j=bestActions.size()-1; j>=0; j--){ //from back to front so we can remove without problems
							if(currentUpperBound - getUpperBoundActionValue(s, bestActions.get(j)).getUpperBound() > epsilon){
								bestActions.remove(j);
							}
						}
					}
				}
			} else {
				if (currentUpperBound<upperBoundValue){
					upperBoundValue = currentUpperBound; //we always remember the minimal upper bound possible
				}
				if(lowerBoundValue - currentLowerBound > epsilon) { //the current action has a value that is smaller by more than epsilon.
					bestActions.clear();
					bestActions.add(i);
					lowerBoundValue = currentLowerBound;
				}
				else if (Math.abs(currentLowerBound - lowerBoundValue)<epsilon){//the current action is epsilon close to the best one
					//another best Action, add it to the set
					bestActions.add(i);
					//We want to stay epsilon close to the best action, so we need to remember the lowest Lower Bound in bestActions, and we need to check that nothing forbidden is in bestActions
					if(lowerBoundValue > currentLowerBound){
						lowerBoundValue = currentLowerBound;
						//The new action is smaller than the old; therefor the epsilon area of the best action shifts, and we have to check, whether we have to remove something from bestActions
						for (int j=bestActions.size()-1; j>=0; j--){ //from back to front so we can remove without problems
							if(getLowerBoundActionValue(s, bestActions.get(j)).getLowerBound() - currentLowerBound > epsilon){
								bestActions.remove(j);
							}
						}
					}

				}
			}
		}
		setQValue(s,new StateValue(lowerBoundValue, upperBoundValue));

		if(bestActions.isEmpty()){System.err.println("Empty set of best actions in StateUpdateSMG. This can never happen.");System.exit(-1);}

		return bestActions;
	}
}