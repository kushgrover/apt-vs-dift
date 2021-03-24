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

import heuristics.HeuristicsSMGModelChecker;
import heuristics.nextstate.HeuristicNextState;
import heuristics.update.StateUpdate;

import explicit.ModelExplorer;
import parser.State;
import prism.PrismException;

public abstract class Heuristic {
	
	protected HeuristicsSMGModelChecker mc;
	protected ModelExplorer pme;
	protected State initialState;
	protected int bound = -1;
	protected HeuristicNextState nextState;
	protected StateUpdate stateUpdate;
	protected boolean verbose;
	protected boolean min = false;
	
	//Statistics
	protected int trialSteps = 0;
	protected int trials = 0;
	
	public Heuristic(HeuristicsSMGModelChecker mc,StateUpdate su, HeuristicNextState nextState, ModelExplorer pme, boolean min) throws PrismException{
		this.mc = mc;
		this.pme = pme;
		this.stateUpdate = su;
		this.nextState = nextState;
		initialState = pme.getDefaultInitialState();
		this.min = min;
	}
	
	public abstract void heuristicStep(State s) throws PrismException;
	public abstract boolean isDone() throws PrismException;
	
	public void compute() throws PrismException{
		State initial = pme.getDefaultInitialState();
		double lowerBoundInitialState = 0;
		double upperBoundInitialState = 1;
		boolean done = false;
		
		stateUpdate.setQValue(initial, new StateValue(lowerBoundInitialState, upperBoundInitialState), 0);
		
		heuristicStart();
		
		while(!done) {
			State s = initial;
			heuristicStep(s);
			done = isDone();
			//mc.getLog().println("Current bounds: " + stateUpdate.getQValue(s).getLowerBound() + " - " + stateUpdate.getQValue(s).getUpperBound());
		}
		reportProgress(trials, trialSteps);
		heuristicStop();
	}
	
	protected void heuristicStart() throws PrismException{
		
	}
	
	protected void heuristicStop() throws PrismException{
		
	}
	
	public int getBound() {
		return bound;
	}
	
	public void setBound(int b) {
		this.bound = b;
	}
	
	public void setVerbose(boolean v) {
		this.verbose = v;
	}
	
	public boolean getVerbose() {
		return this.verbose;
	}
	
	public StateValue getInitialStateValue() {
		return stateUpdate.getQValue(initialState,0) == null ? new StateValue(-1, -1) : stateUpdate.getQValue(initialState,0);
	}
	
	protected void reportProgress(int trials, int trialSteps) {
		if(verbose) {
			double lowerBoundInitialState = stateUpdate.getQValue(initialState,0).getLowerBound();
			double upperBoundInitialState = stateUpdate.getQValue(initialState,0).getUpperBound();
			mc.getLog().println("Trials: " + trials);
			mc.getLog().println("Trial steps: " + trialSteps);
			mc.getLog().println("Trial steps avg: " + (double)trialSteps/(double)trials);
			mc.getLog().println("Explore: " + stateUpdate.getNumExplored());
			mc.getLog().println("Lower bound: " + lowerBoundInitialState);
			mc.getLog().println("Upper bound: " + upperBoundInitialState);
			//mc.getLog().println("All values: " + stateUpdate.getNumAllValues());
			mc.getLog().flush();
		}
	}
}
