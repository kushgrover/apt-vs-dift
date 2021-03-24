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
import java.util.List;

import heuristics.HeuristicsSMGModelChecker;
import heuristics.nextstate.HeuristicNextState;
import heuristics.nextstate.HeuristicNextStateFactory;
import heuristics.update.StateUpdate;
import heuristics.update.StateValueContainer;

import explicit.ModelExplorer;

import parser.State;
import prism.PrismException;

public abstract class HeuristicParallel extends Heuristic{

	private List<Thread> workers = new ArrayList<Thread>();
	
	public HeuristicParallel(HeuristicsSMGModelChecker mc, StateUpdate su, HeuristicNextState nextState, ModelExplorer pme, boolean min) throws PrismException{
		super(mc, su, nextState, pme, min);
	}
	
	@Override
	protected void heuristicStart() throws PrismException{
		int threads = Runtime.getRuntime().availableProcessors()-1;
		if(threads == 0) threads = 1;
		startThreads(threads);
	}
	
	@Override
	protected void heuristicStop() throws PrismException{
		threadsInterrupt();
	}

	private void startThreads(int threads) throws PrismException{
		mc.getLog().println("Starting threads: " + threads);
		for(int i=0;i<threads;i++) {
			Thread worker = createWorker(stateUpdate, nextState, initialState);
			workers.add(worker);
			worker.start();
		}
	}
	
	private void threadsInterrupt() {
		for(int i=0;i<workers.size();i++) {
			workers.get(i).interrupt();
		}
	}
	
	protected abstract HeuristicWorker createWorker(StateUpdate suWorker, HeuristicNextState nsWorker, parser.State isWorkers) throws PrismException;
	
	protected static abstract class HeuristicWorker extends Thread {

		protected StateUpdate stateUpdateWorker;
		protected parser.State initialStateWorker;
		protected HeuristicNextState nextStateWorker;
		protected ModelExplorer modelExplorerWorker;
		
		public HeuristicWorker(StateUpdate suWorker, HeuristicNextState nsWorker, parser.State isWorker) throws PrismException{
			this.initialStateWorker = isWorker;
			StateValueContainer svc = suWorker.getStateValueContainer();
			//this.modelExplorerWorker = suWorker.getModelExplorer().clone(); wont work anymore, cause I disabled clone to exist, cause it didnt exist in the cachedModelExplorer
			this.stateUpdateWorker = suWorker.clone();
			this.stateUpdateWorker.setStateValueContainer(svc);
			try {
				this.nextStateWorker = HeuristicNextStateFactory.getHeuristicNextState(nsWorker.getType(), this.modelExplorerWorker, this.stateUpdateWorker);
				this.stateUpdateWorker.setModelExplorer(this.modelExplorerWorker);
				this.stateUpdateWorker.setCoalition(suWorker.getCoalition());
				this.stateUpdateWorker.setTarget(suWorker.getTarget());
			} catch (PrismException e) {
				System.out.println(e);
				e.printStackTrace();
				System.exit(1);
			}
		}
		
		public void run() {
			while(!Thread.currentThread().isInterrupted()) {
				try {
					heuristicStep();
				} catch (PrismException e) {
					System.out.println(e);
					e.printStackTrace();
					System.exit(1);
				}
			}
		}
		
		public abstract void heuristicStep() throws PrismException;
	}
}