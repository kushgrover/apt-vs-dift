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

public class HeuristicRTDPParallel extends HeuristicParallel{

	private HeuristicRTDP heuristicRTDP;
	
	public HeuristicRTDPParallel(HeuristicsSMGModelChecker mc, StateUpdate su, HeuristicNextState nextState, ModelExplorer pme, boolean min) throws PrismException{
		super(mc, su, nextState, pme, min);
		this.heuristicRTDP = new HeuristicRTDP(mc, su, nextState, pme, min);
		this.heuristicRTDP.setVerbose(getVerbose());
		this.heuristicRTDP.setBound(getBound());
	}
	
	@Override
	protected HeuristicWorker createWorker(StateUpdate suWorker, HeuristicNextState nsWorker, State isWorkers) throws PrismException{
		return new RTDPWorker(suWorker, nsWorker, isWorkers);
	}

	@Override
	public void heuristicStep(State s) throws PrismException {
	
	}

	public boolean isDone() throws PrismException {
		return this.heuristicRTDP.isDone();
	}
	
	protected class RTDPWorker extends HeuristicWorker {

		private HeuristicRTDP heuristicRTDPWorker; 
		
		public RTDPWorker(StateUpdate suWorker, HeuristicNextState nsWorker, parser.State isWorker) throws PrismException{
			super(suWorker, nsWorker, isWorker);
			this.heuristicRTDPWorker = new HeuristicRTDP(mc, stateUpdateWorker, nextStateWorker, modelExplorerWorker, min);
			this.heuristicRTDPWorker.setVerbose(getVerbose());
			this.heuristicRTDPWorker.setBound(getBound());
		}

		@Override
		public void heuristicStep() throws PrismException {
			this.heuristicRTDPWorker.heuristicStep(initialStateWorker);
		}
	}
}