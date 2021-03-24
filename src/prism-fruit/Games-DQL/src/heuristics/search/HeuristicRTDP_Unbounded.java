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

package heuristics.search;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Stack;
import java.util.concurrent.ExecutorCompletionService;

import heuristics.CachedModelExplorer;
import heuristics.HeuristicsSMGModelChecker;
import heuristics.nextstate.HeuristicNextState;
import heuristics.search.HeuristicLRTDP.VisitedState;
import heuristics.update.StateUpdate;

import explicit.ConstructModel;
import explicit.Distribution;
import explicit.MDP;
import explicit.MDPModelChecker;
import explicit.MDPSimple;
import explicit.Model;
import explicit.ModelExplorer;
import explicit.StateValues;

import parser.State;
import parser.ast.Expression;
import prism.ECComputer;
import prism.PrismException;

/**
 * Only kept for reference, outdated, unused; Maxi, 12.01.18
 */
public class HeuristicRTDP_Unbounded {//extends Heuristic{

//	private long modelCheckingTime = 0;
//	private boolean precomp = false;
//
//	public HeuristicRTDP_Unbounded(HeuristicsSMGModelChecker mc, StateUpdate su, HeuristicNextState nextState, ModelExplorer pme, boolean min) throws PrismException{
//		super(mc, su, nextState, pme, min);
//		if(precomp) {
//			preComp();
//		}
//	}
//
//	@Override
//	protected void heuristicStart() throws PrismException{
//		modelCheckingTime = System.currentTimeMillis();
//	}
//
//	@Override
//	protected void heuristicStop() throws PrismException{
//		long duration = System.currentTimeMillis() - modelCheckingTime;
//		mc.getLog().println();
//		mc.getLog().println("Heuristic model checking time in " + ((double)duration/1000)  + " secs.");
//	}
//
//	private Set<State> seen = new HashSet<State>();
//	private int buildPartialEvery = 100000;//default 100000
//	private int precompCount = 0;
//	private int precompThreshold = 15;
//
//	@Override
//	public void heuristicStep(State s) throws PrismException {
//		Stack<State> visited = new Stack<State>();
//		while(!stateUpdate.isTarget(s) && !stateUpdate.isZero(s) && !visited.contains(s)) {
//			seen.add(s);
//			visited.add(s);
//			int bestAction = stateUpdate.update(s);
//			s = nextState.sample(s,bestAction);
//			if(!precomp && trialSteps % buildPartialEvery == 0) {
//				updatePrecomp();
//				precompCount++;
//				if(precompCount > precompThreshold) {
//					buildPartialEvery = buildPartialEvery*100;
//					precompCount = 0;
//				}
//			}
//			trialSteps++;
//		}
//		seen.add(s);
//		visited.add(s);
//		stateUpdate.setTargetOrZeroValue(s);
//
//		while(!visited.empty()) {
//			State vs = visited.pop();
//			stateUpdate.update(vs);
//		}
//
//		if(trials % 10000 == 0) {
//			reportProgress(trials, trialSteps);
//		}
//		trials++;
//	}
//
//	@Override
//	public boolean isDone() throws PrismException {
//		StateValue sv = stateUpdate.getQValue(initialState);
//		if(sv != null) {
//			double lowerBoundInitialState = sv.getLowerBound();
//			double upperBoundInitialState = sv.getUpperBound();
//			return stateUpdate.same(lowerBoundInitialState, upperBoundInitialState);
//		}
//		return false;
//	}
//
//	private void updatePrecomp() throws PrismException{
//		MDPModelChecker mdpModelChecker = getMC();
//		MDPSimple mdp = buildPartialModel();
//		mdp.findDeadlocks(true);
//		BitSet target = computeTarget(mdpModelChecker, mdp);
//
//		if(min) {
//			BitSet t = new BitSet();
//			t.or(target);
//			t.or(unvisited);
//			BitSet prob0 = computeProb0(mdpModelChecker, mdp, t);
//			for (int i = prob0.nextSetBit(0); i >= 0; i = prob0.nextSetBit(i+1)) {
//				State s = mdp.getStatesList().get(i);
//				stateUpdate.setQValue(s, new StateValue(0, 0));
//				stateUpdate.setZero(s, true);
//			}
//			BitSet prob1 = computeProb1(mdpModelChecker, mdp, target);
//			for (int i = prob1.nextSetBit(0); i >= 0; i = prob1.nextSetBit(i+1)) {
//				State s = mdp.getStatesList().get(i);
//				if(seen.contains(s)) {
//					stateUpdate.setQValue(s, new StateValue(1, 1));
//					stateUpdate.setTarget(s, true);
//				}
//			}
//		} else {
//			target.or(unvisited);
//			BitSet prob0 = computeProb0(mdpModelChecker, mdp, target);
//			for (int i = prob0.nextSetBit(0); i >= 0; i = prob0.nextSetBit(i+1)) {
//				State s = mdp.getStatesList().get(i);
//				if(seen.contains(s)) {
//					stateUpdate.setQValue(s, new StateValue(0, 0));
//					stateUpdate.setZero(s, true);
//				}
//			}
//			/*BitSet prob1 = computeProb1(mdpModelChecker, mdp, target);
//			for (int i = prob1.nextSetBit(0); i >= 0; i = prob1.nextSetBit(i+1)) {
//				State s = mdp.getStatesList().get(i);
//				if(seen.contains(s)) {
//					stateUpdate.setQValue(s, new StateValue(1, 1));
//					stateUpdate.setTarget(s, true);
//				}
//			}*/
//			collapseMECs(mdp);
//		}
//	}
//
//	private MDPSimple mdp = new MDPSimple();
//	private Map<State, Integer> state2Index = new HashMap<State, Integer>();
//	private BitSet unvisited = new BitSet();
//
//	private MDPSimple buildPartialModel() throws PrismException{
//		List<State> statesList = mdp.getStatesList();
//		if(statesList == null) {
//			statesList = new ArrayList<State>();
//		}
//		unvisited.clear();
//		if(mdp.getFirstInitialState() == -1) {
//			int index = mdp.addState();
//			mdp.addInitialState(index);
//			state2Index.put(initialState, index);
//			statesList.add(initialState);
//		}
//		Iterator<State> it = seen.iterator();
//		while(it.hasNext()) {
//			State s = it.next();
//			Integer index = state2Index.get(s);
//			if(index == null) {
//				index = mdp.addState();
//				state2Index.put(s, index);
//				statesList.add(s);
//			}
//			List<Distribution> dists = buildAllDistributionsIn(mdp, s, seen, statesList);
//			mdp.clearState(index);
//			for(int i=0;i<dists.size();i++) {
//				mdp.addChoice(index, dists.get(i));
//			}
//		}
//		mdp.setStatesList(statesList);
//		mc.getLog().println("Model built " + mdp.getNumStates());
//		mc.getLog().flush();
//		return mdp;
//	}
//
//	private List<Distribution> buildAllDistributionsIn(MDPSimple mdp, State s, Set<State> seen, List<State> statesList) throws PrismException{
//		List<Distribution> dists = new ArrayList<Distribution>();
//		pme.queryState(s);
//		int choices = pme.getNumChoices();
//		for(int i=0;i<choices;i++) {
//			Distribution d = new Distribution();
//			int trans = pme.getNumTransitions(i);
//			boolean inSeen = false;
//			for(int j=0;j<trans;j++) {
//				double prob = pme.getTransitionProbability(i, j);
//				State t = pme.computeTransitionTarget(i,j);
//				inSeen = inSeen || seen.contains(t);
//				Integer index = state2Index.get(t);
//				if(index == null) {
//					index = mdp.addState();
//					state2Index.put(t, index);
//					statesList.add(t);
//				}
//				d.add(index, prob);
//				if(!seen.contains(t)) {
//					unvisited.set(index);
//				}
//			}
//			dists.add(d);
//		}
//		return dists;
//	}
//
//
//	private void preComp() {
//		try {
//			ConstructModel constructModel = new ConstructModel(mc, mc.getEngine(), null, null);
//			MDP mdp = (MDP)constructModel.constructModel(mc.getModulesFile());
//			MDPModelChecker mdpModelChecker = getMC();
//			BitSet target = computeTarget(mdpModelChecker, mdp);
//			BitSet zero = computeProb0(mdpModelChecker, mdp, target);
//			BitSet one = computeProb1(mdpModelChecker, mdp, target);
//
//			List<State> statesList = mdp.getStatesList();
//			for(int i=0;i<statesList.size();i++) {
//				State s = statesList.get(i);
//				if(target.get(i) || one.get(i)) {
//					stateUpdate.setTarget(s, true);
//				} else {
//					stateUpdate.setTarget(s, false);
//				}
//				if(zero.get(i)) {
//					stateUpdate.setZero(s, true);
//				}
//				pme.queryState(s);
//			}
//			if(!min) {
//				collapseMECs(mdp);
//			}
//		} catch(PrismException e) {
//			mc.getLog().print(e);
//		}
//	}
//
//	private BitSet computeProb0(MDPModelChecker mdpModelChecker, MDP mdp, BitSet target) throws PrismException{
//		return mdpModelChecker.prob0(mdp, null, target, min, null);
//	}
//
//	private BitSet computeProb1(MDPModelChecker mdpModelChecker, MDP mdp, BitSet target) throws PrismException{
//		return mdpModelChecker.prob1(mdp, null, target, min, null);
//	}
//
//	private BitSet computeTarget(MDPModelChecker mdpModelChecker, MDP mdp) throws PrismException {
//		Expression targetExp = null;
//		if(mc.getExpression().getOperand1() != null) {
//			targetExp = mc.getExpression();
//		} else {
//			targetExp = mc.getExpression().getOperand2();
//		}
//		return mdpModelChecker.checkExpression(mdp, targetExp).getBitSet();
//	}
//
//	private MDPModelChecker getMC() throws PrismException{
//		MDPModelChecker mdpModelChecker = new MDPModelChecker(mc);
//		mdpModelChecker.setModulesFileAndPropertiesFile(mc.getModulesFile(), mc.getPropertiesFile());
//		return mdpModelChecker;
//	}
//
//	private void collapseMECs(MDP mdp) throws PrismException {
//		mc.getLog().println("Starting MEC collapsing");
//		mc.getLog().flush();
//		long start = System.currentTimeMillis();
//		BitSet all = new BitSet();
//		BitSet notOne = new BitSet();
//		for(int i=0;i<mdp.getNumStates();i++) {
//			all.set(i);
//			StateValue sv = stateUpdate.getQValue(mdp.getStatesList().get(i));
//			if(sv != null) {
//				if(sv.getUpperBound() < 1.0) {
//					notOne.set(i);
//				}
//			}
//		}
//		all.xor(notOne);
//		explicit.ECComputer ec = explicit.ECComputer.createECComputer(mc, mdp);
//		ec.computeMECStates(all);
//		List<BitSet> mecs = ec.getMECStates();
//		for(int i=0;i<mecs.size();i++) {
//			BitSet mec = mecs.get(i);
//			List<Map<State, Double>> actions = getAllLeavingMEC(mdp, mec);
//			if(actions.size() > 0) {
//				if(!min && containsTarget(mdp, mec)) {
//					for (int j = mec.nextSetBit(0); j >= 0; j = mec.nextSetBit(j+1)) {
//						State s = mdp.getStatesList().get(j);
//						stateUpdate.setQValue(s, new StateValue(1.0, 1.0));
//						stateUpdate.setTarget(s, true);
//					}
//				} else {
//					for (int j = mec.nextSetBit(0); j >= 0; j = mec.nextSetBit(j+1)) {
//						collapse(mdp, j, actions);
//					}
//				}
//			} else {
//				if(!min && containsTarget(mdp, mec)) {
//					for (int j = mec.nextSetBit(0); j >= 0; j = mec.nextSetBit(j+1)) {
//						State s = mdp.getStatesList().get(j);
//						stateUpdate.setQValue(s, new StateValue(1, 1));
//						stateUpdate.setTarget(s, true);
//					}
//				}
//			}
//		}
//		long duration = System.currentTimeMillis() - start;
//		mc.getLog().println("MEC collapsing done " + (double)duration/1000 + " secs.");
//		mc.getLog().flush();
//	}
//
//	private boolean containsTarget(MDP mdp, BitSet mec) throws PrismException{
//		for (int i = mec.nextSetBit(0); i >= 0; i = mec.nextSetBit(i+1)) {
//			State s = mdp.getStatesList().get(i);
//			if(stateUpdate.isTarget(s)) {
//				return true;
//			}
//		}
//		return false;
//	}
//
//	private void collapse(MDP mdp, int s, List<Map<State, Double>> actions) throws PrismException{
//		CachedModelExplorer cme = (CachedModelExplorer)pme;
//		cme.updateCache(mdp.getStatesList().get(s), actions);
//	}
//
//	private List<Map<State, Double>> getAllLeavingMEC(MDP mdp, BitSet mec) {
//		List<Map<State, Double>> actions = new ArrayList<Map<State, Double>>();
//		for (int s = mec.nextSetBit(0); s >= 0; s = mec.nextSetBit(s+1)) {
//			for(int i = 0; i < mdp.getNumChoices(s);i++) {
//				boolean all = mdp.allSuccessorsInSet(s, i, mec);
//				if(!all) {
//					actions.add(getDistribution(mdp, s, i));
//				}
//			}
//		}
//		return actions;
//	}
//
//	private Map<State, Double> getDistribution(MDP mdp, int s, int i) {
//		Map<State, Double> d = new HashMap<State, Double>();
//		Iterator<Entry<Integer,Double>> it =  mdp.getTransitionsIterator(s, i);
//		while(it.hasNext()) {
//			Entry<Integer, Double> e = it.next();
//			d.put(mdp.getStatesList().get(e.getKey()), e.getValue());
//		}
//		return d;
//	}
}