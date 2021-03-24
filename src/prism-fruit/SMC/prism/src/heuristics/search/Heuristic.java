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
import parser.State;
import prism.ModelGenerator;
import prism.PrismException;

public abstract class Heuristic
{
	protected final boolean min;
	protected HeuristicsSMGModelChecker smgModelChecker;
	protected ModelGenerator modelGenerator;
	protected int rewardIndex;
	protected State initialState;
	protected int bound = -1;
	protected HeuristicNextState nextState;
	protected StateUpdate stateUpdate;
	protected boolean verbose;
	protected double epsilon;

	//Statistics
	protected int trialSteps = 0;
	protected int trials = 0;

	public Heuristic(HeuristicsSMGModelChecker smgModelChecker, StateUpdate stateUpdate, HeuristicNextState nextState, ModelGenerator modelGenerator,
			boolean min)
			throws PrismException
	{
		this.smgModelChecker = smgModelChecker;
		this.modelGenerator = modelGenerator;
		this.stateUpdate = stateUpdate;
		this.nextState = nextState;
		this.initialState = modelGenerator.getInitialState();
		this.min = min;
	}

	public void setRewardIndex(int rewardIndex)
	{
		this.rewardIndex = rewardIndex;
	}

	public abstract void exploreAndUpdate(State s) throws PrismException;

	public abstract boolean isDone() throws PrismException;

	/**
	 * The learning algorithm for computing max/min prob
	 *
	 * @throws PrismException
	 */
	public void computeProb() throws PrismException
	{
		State initial = modelGenerator.getInitialState();

		if (stateUpdate.getBound() == -1) {
			stateUpdate.setQValue(initial, new StateValue(0d, 1d));
		}
		else {
			stateUpdate.setQValue(initial, new StateValue(0d, 1d), 0);
		}

		heuristicStart();
		do {
			exploreAndUpdate(initial);
		} while (!isDone());
 		reportProgress(trials, trialSteps);
		heuristicStop();
	}

	protected abstract void heuristicStart() throws PrismException;

	protected abstract void heuristicStop() throws PrismException;

	public int getBound()
	{
		return bound;
	}

	public void setBound(int bound)
	{
		this.bound = bound;
	}

	public void setVerbose(boolean v)
	{
		this.verbose = v;
	}

	public boolean getVerbose()
	{
		return this.verbose;
	}

	public void setEpsilon(double epsilon)
	{
		this.epsilon = epsilon;
	}

	public StateValue getInitialStateValue()
	{
		final StateValue qValue = stateUpdate.getQValue(initialState);
		if (qValue == null) {
			return new StateValue(-1.0, -1.0);
		}
		return qValue;
	}

	protected void reportProgress(int trials, int trialSteps)
	{
		int depth = getBound() == -1 ? -1 : 0;
		if (verbose) {
			double lowerBoundInitialState = stateUpdate.getQValue(initialState, depth).lowerBound;
			double upperBoundInitialState = stateUpdate.getQValue(initialState, depth).upperBound;
			smgModelChecker.getLog().println("Trials: " + trials);
			smgModelChecker.getLog().println("Trial steps: " + trialSteps);
			smgModelChecker.getLog().println("Trial steps avg: " + (double) trialSteps / (double) trials);
			smgModelChecker.getLog().println("Explore: " + stateUpdate.getNumExplored());
			smgModelChecker.getLog().println("Lower bound: " + lowerBoundInitialState);
			smgModelChecker.getLog().println("Upper bound: " + upperBoundInitialState);
			//smgModelChecker.getLog().println("All values: " + stateUpdate.getNumAllValues());
			smgModelChecker.getLog().flush();
		}
	}
}
