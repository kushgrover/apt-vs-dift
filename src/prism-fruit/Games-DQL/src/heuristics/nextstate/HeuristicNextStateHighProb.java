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

import java.util.Random;

import explicit.ModelExplorer;

import parser.State;
import prism.PrismException;


public class HeuristicNextStateHighProb extends HeuristicNextState{

	public HeuristicNextStateHighProb(ModelExplorer pme, StateUpdate stateUpdate) throws PrismException{
		super(pme, stateUpdate);
	}
	
	public State sample(State s, int action, int depth) throws PrismException{
		pme.queryState(s);
		int trans = pme.getNumTransitions(action);
		double probs[] = new double[trans];
		for(int i=0;i<trans;i++) {
			probs[i] = pme.getTransitionProbability(action, i);
		}
		return sampleFromDist(action, probs);
		
	}
	
	public State sample(State s, int action) throws PrismException{
		return sample(s, action, -1);
	}
	
	public NextState getType() {
		return NextState.HIGH_PROB;
	}
				
}
