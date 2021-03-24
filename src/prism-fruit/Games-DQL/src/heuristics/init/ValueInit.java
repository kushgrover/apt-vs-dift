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

package heuristics.init;

import heuristics.update.StateUpdate;
import explicit.ModelExplorer;
import parser.State;
import prism.PrismException;

public class ValueInit {
	
	private StateUpdate stateUpdate = null;
	
	public ValueInit(StateUpdate su) {
		this.stateUpdate = su;
	}
	
	//Add self loop detection

	public double getLowerBoundValue(ModelExplorer pme, State s, int depth) throws PrismException{
		if(depth == -1) {
			if(stateUpdate.isTarget(s)) return 1;
			if(stateUpdate.isZero(s)) return 0;
			return 0;
		}
		checkDepth(depth);
		
		if(stateUpdate.isTarget(s)) return 1;
		if(depth >= stateUpdate.getBound()) return 0;
		
		return 0;
	}
	
	public double getUpperBoundValue(ModelExplorer pme, State s, int depth) throws PrismException{
		if(depth == -1) {
			if(stateUpdate.isTarget(s)) return 1;
			if(stateUpdate.isZero(s)) return 0;
			return 1;
		}
		checkDepth(depth);
		
		if(stateUpdate.isTarget(s)) return 1;
		if(depth >= stateUpdate.getBound()) return 0;
		if(stateUpdate.isSelfLoop(s)) return 0;
		
		return 1;
	}
	
	
	
	private void checkDepth(int depth) throws PrismException {
		if(depth > stateUpdate.getBound()+1) {
			throw new PrismException("Depth cannot be bigger than bound");
		}
	}

}
