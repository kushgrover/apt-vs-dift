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
import heuristics.search.HeuristicParallel.HeuristicWorker;
import heuristics.search.HeuristicRTDPParallel.RTDPWorker;
import heuristics.update.StateUpdate;

import java.util.BitSet;
import java.util.HashMap;
import java.util.Map;
import java.util.Stack;

import explicit.ModelExplorer;

import parser.State;
import prism.PrismException;

//Labelled version of LRTDP

public class HeuristicLRTDPParallel extends HeuristicParallel{

	private HeuristicLRTDP heuristicLRTDP;
	
	public HeuristicLRTDPParallel(HeuristicsSMGModelChecker mc, StateUpdate su, HeuristicNextState nextState, ModelExplorer pme, boolean min) throws PrismException{
		super(mc, su, nextState, pme, min);
		this.heuristicLRTDP = new HeuristicLRTDP(mc, su, nextState, pme, min);
		this.heuristicLRTDP.setVerbose(getVerbose());
		this.heuristicLRTDP.setBound(getBound());
	}

	@Override
	protected HeuristicWorker createWorker(StateUpdate suWorker, HeuristicNextState nsWorker, State isWorkers) throws PrismException{
		return new LRTDPWorker(suWorker, nsWorker, isWorkers);
	}

	@Override
	public void heuristicStep(State s) throws PrismException {
		
	}

	@Override
	public boolean isDone() throws PrismException {
		double lowerBound = stateUpdate.getQValue(initialState,0).getLowerBound();
		double upperBound = stateUpdate.getQValue(initialState,0).getUpperBound();
		return stateUpdate.same(lowerBound, upperBound);
	}
	
	protected class LRTDPWorker extends HeuristicWorker {

		private HeuristicLRTDP heuristicLRTDPWorker; 
		
		public LRTDPWorker(StateUpdate suWorker, HeuristicNextState nsWorker, parser.State isWorker) throws PrismException{
			super(suWorker, nsWorker, isWorker);
			this.heuristicLRTDPWorker = new HeuristicLRTDP(mc, stateUpdateWorker, nextStateWorker, modelExplorerWorker, min);
			this.heuristicLRTDPWorker.setVerbose(getVerbose());
			this.heuristicLRTDPWorker.setBound(getBound());
		}

		@Override
		public void heuristicStep() throws PrismException {
			this.heuristicLRTDPWorker.heuristicStep(initialStateWorker);
		}
	}
}
