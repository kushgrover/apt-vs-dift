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
import java.util.Map.Entry;

import explicit.ModelExplorer;
import explicit.SMG;
import explicit.STPGExplicit;

import parser.State;
import parser.Values;
import parser.ast.Expression;
import parser.ast.ExpressionIdent;
import parser.ast.LabelList;
import parser.ast.ModulesFile;
import prism.ModelType;
import prism.Pair;
import prism.PrismException;
import prism.PrismUtils;
import heuristics.init.ValueInit;
import heuristics.search.StateValue;

public abstract class StateUpdate implements Cloneable {

	private static final Random random = new Random();
	protected ModelExplorer pme;
	private Values constantValues;
	private Expression target;
	private State initialState;
	private StateValueContainer stateValueContainer;
	private ValueInit valueInit = new ValueInit(this);
	protected int bound = -1;
	private LabelList specialLabels;
	private Set<Integer> coalition;
	protected boolean min;
	protected double epsilon;

	public StateUpdate(ModelExplorer me, StateValueContainer container, int bound, boolean min, double epsilon) throws PrismException{
		this.pme = me;
		this.bound = bound;
		this.specialLabels = new LabelList();
		this.specialLabels.addLabel(new ExpressionIdent("deadlock"), new ExpressionIdent("deadlock"));
		this.specialLabels.addLabel(new ExpressionIdent("init"), new ExpressionIdent("init"));
		this.stateValueContainer = container;
		if(me != null) {
			this.initialState = pme.getDefaultInitialState();
		}
		this.min = min;
		this.epsilon = epsilon;
	}
	
	public abstract int update(State s, int depth) throws PrismException;
	
	public abstract int update(State s, boolean isMaxState) throws PrismException;

	//Relic of an optimization that did not work well. Maybe look at this again in the future. games-brtdp gitlog should contain the necessary information
	//public abstract int update(State s, HashMap<State,BitSet> state2leavingActions) throws PrismException;


	public int sampleBestAction(State s, boolean isMaxState) throws PrismException{
		ArrayList<Integer> bestActions = getBestActions(s, isMaxState).first;
		//pick a random action from the set of best actions
		int bestAction = bestActions.get(random.nextInt(bestActions.size()));
		return bestAction;
	}

	public abstract Pair<ArrayList<Integer>,StateValue> getBestActions(State s, boolean isMaxState) throws PrismException;
	
	private Map<State, Boolean> cacheState2Target = new HashMap<State,Boolean>();
	private Map<State, Boolean> cacheState2SelfLoop = new HashMap<State,Boolean>();
	private Map<State, Boolean> cacheState2Zero = new HashMap<State,Boolean>();
	
	public boolean isTarget(State s) throws PrismException{
		Boolean cachedValue = cacheState2Target.get(s) ;
		if(cachedValue == null) {
			pme.queryState(s);
			Expression evTarget = target.deepCopy();
			specialLabels.setLabel(0, pme.getNumTransitions() == 0 ? Expression.True() : Expression.False());
			specialLabels.setLabel(1, initialState.equals(s) ? Expression.True() : Expression.False());
			cachedValue = evTarget.evaluateBoolean(constantValues, s);
			cacheState2Target.put(s, cachedValue);
		} 
		return cachedValue;
	}
	
	public boolean isSelfLoop(State s) throws PrismException{
		Boolean cachedValue = cacheState2SelfLoop.get(s);
		if(cachedValue == null) {
			pme.queryState(s);
			int numChoices = pme.getNumChoices();
			for(int i=0;i<numChoices;i++) {
				int trans = pme.getNumTransitions(i);
				for(int j=0;j<trans;j++) {
					State succ = pme.computeTransitionTarget(i,j);
					if(!succ.equals(s)) {
						cacheState2SelfLoop.put(s,false);
						return false;
					}
				}
			}
			cachedValue = true;
			cacheState2SelfLoop.put(s, cachedValue);
		}
		return cachedValue;
	}
	
	public boolean isZero(State s) throws PrismException{
		Boolean cachedValue = cacheState2Zero.get(s);
		if(cachedValue == null) {
			return false;
		}
		return true;
	}
	
	public boolean visited(State s) {
		return getQValue(s) == null ? true : false;
	}
	
	public void setTarget(State s, boolean t) {
		cacheState2Target.put(s, t);
		if(t) {
			setQValue(s, new StateValue(1, 1));
		}
	}
	
	public void setZero(State s, boolean t) {
		cacheState2Zero.put(s, t);
		if(t) {
			setQValue(s, new StateValue(0, 0));
		}
	}
	
	public int getBound() {
		return bound;
	}
	
	public void setQValue(State s, StateValue val, int depth) {
		stateValueContainer.setQValue(s, val, depth);
	}
	
	public StateValue getQValue(State s, int depth) {
		return stateValueContainer.getQValue(s, depth);
	}
	
	public void setQValue(State s, StateValue val) {
		stateValueContainer.setQValue(s, val, 0);
	}
	
	public StateValue getQValue(State s) {
		return stateValueContainer.getQValue(s, 0);
	}
	
	public int getCurrentAction(State s) {
		Integer cA = stateValueContainer.getCurrentAction(s);
		if(cA == null) {
			return 0;
		} 
		return cA;
	}
	
	public void setCurrentAction(State s, int a) {
		stateValueContainer.setCurrentAction(s, a);
	}
	
	public void setTargetOrBoundValue(State s, int depth) throws PrismException{
		if(isTarget(s)) {
			setQValue(s, new StateValue(1.0, 1.0), depth);
			return;
		}
		if(depth >= bound) {
			setQValue(s, new StateValue(0, 0), depth);
			return;
		}
	}
	
	public void setTargetOrZeroValue(State s) throws PrismException{
		if(isTarget(s)) {
			setQValue(s, new StateValue(1.0, 1.0));
			return;
		}
		if(isZero(s)) {
			setQValue(s, new StateValue(0, 0));
			return;
		}
	}
	
	public int getStateValueNumber() {
		int stateValuesNumber = 0;
		/*Iterator<State> it = qValues.keySet().iterator();
		while(it.hasNext()) {
			State s = it.next();
			stateValuesNumber += qValues.get(s).size();
			if(qValues.get(s).size() == 0) {
				System.out.println("s " + s);
				System.exit(1);
			}
		}*/
		return stateValuesNumber;
	}
	
	public int getNumExplored() {
		return stateValueContainer.getSize();
	}
	
	public int getNumAllValues() {
		return stateValueContainer.getNumAllValues();
	}
	
	public boolean same(double d1, double d2) {
		return PrismUtils.doublesAreClose(d1, d2, epsilon, true);
	}
	
	public void setConstantValues(Values cV) {
		constantValues = cV;
	}
	
	public void setTarget(Expression e) {
		this.target = e;
	}
	
	public Expression getTarget() {
		return this.target;
	}
	
	public void setCoalition(Set<Integer> coalition) {
		this.coalition = coalition;
	}
	
	public Set<Integer> getCoalition() {
		return this.coalition;
	}
	
	public int getPlayer(State s) throws PrismException{
		pme.queryState(s);
		return coalition.contains(pme.getPlayerForState()) ? 1 : 2;
	}
	
	public void setModelExplorer(ModelExplorer me) throws PrismException{
		this.pme = me;
		this.initialState = pme.getDefaultInitialState();
	}
	
	public void setStateValueContainer(StateValueContainer container) {
		this.stateValueContainer = container;
	}
	
	public ModelExplorer getModelExplorer() {
		return this.pme;
	}

	public double getEpsilon(){
		return epsilon;
	}

	public StateValueContainer getStateValueContainer() {
		return this.stateValueContainer;
	}
	
	protected StateValue getLowerBoundActionValue(State s, int a, int depth) throws PrismException {
		pme.queryState(s);
		int trans = pme.getNumTransitions(a);
		double bestVal = 0;
		for(int i=0;i<trans;i++) {
			pme.queryState(s);
			double prob = pme.getTransitionProbability(a, i);
			State succ = pme.computeTransitionTarget(a, i);
			StateValue tVal = null;
			if(depth == -1) {
				tVal = getQValue(succ);
			} else {
				tVal = getQValue(succ, depth);
			}
			
			if(tVal != null) {
				bestVal += prob * tVal.getLowerBound();
			} else {
				double succBestValue = valueInit.getLowerBoundValue(pme, succ, depth);
				bestVal += prob * succBestValue;
			}
		}
		return new StateValue(bestVal, -1);
	}
	
	protected StateValue getUpperBoundActionValue(State s, int a, int depth) throws PrismException {
		pme.queryState(s);
		int trans = pme.getNumTransitions(a);
		double worseVal = 0;
		for(int i=0;i<trans;i++) {
			pme.queryState(s);
			double prob = pme.getTransitionProbability(a, i);
			State succ = pme.computeTransitionTarget(a, i);
			StateValue tVal = null;
			if(depth == -1) {
				tVal = getQValue(succ);
			} else {
				tVal = getQValue(succ, depth);
			}
			if(tVal != null) {
				worseVal += prob * tVal.getUpperBound();
			} else {
				double succWorseValue = valueInit.getUpperBoundValue(pme, succ, depth);
				worseVal += prob * succWorseValue;
			}
		}
		return new StateValue(-1, worseVal);
	}
	
	public StateValue getLowerBoundActionValue(State s, int a) throws PrismException {
		return getLowerBoundActionValue(s, a, -1);
	}
	
	public StateValue getUpperBoundActionValue(State s, int a) throws PrismException {
		return getUpperBoundActionValue(s, a, -1);
	}
	
	public abstract StateUpdate clone();
	
}
