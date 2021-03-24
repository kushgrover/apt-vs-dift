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

import java.util.Arrays;
import java.util.Random;

import heuristics.nextstate.HeuristicNextStateFactory.NextState;
import heuristics.update.StateUpdate;
import explicit.ModelExplorer;

import parser.State;
import prism.PrismException;


public abstract class HeuristicNextState implements Cloneable{

	protected ModelExplorer pme;
	protected StateUpdate stateUpdate;
	protected Random random = new Random();
	
	public HeuristicNextState(ModelExplorer pme, StateUpdate stateUpdate) throws PrismException{
		this.pme = pme;
		this.stateUpdate = stateUpdate;
	}
	
	public abstract State sample(State s, int action, int depth) throws PrismException;
	
	public abstract State sample(State s, int action) throws PrismException;
	
	public abstract NextState getType();
	
	protected State sampleFromDist(int action, double[] dist) throws PrismException {
		for(int i=1;i<dist.length;i++) {
			dist[i] =  dist[i-1] + dist[i];
		}

		final double randomValue = random.nextDouble();

		int index = Arrays.binarySearch(dist, randomValue);

		if (index < 0) {
			index = -index - 1;
		}

		if (index >= 0 &&
				index < dist.length &&
				randomValue < dist[index]) {
			return pme.computeTransitionTarget(action, index);
		}

		//throw new PrismException("When would this happen? Dist unset in sampleFromDist or sth like that. Happens for white, apparently.");
		return pme.computeTransitionTarget(action, random.nextInt(dist.length));
	}
	
}
