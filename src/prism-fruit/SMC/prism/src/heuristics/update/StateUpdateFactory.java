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

import heuristics.CachedModelGenerator;
import prism.ModelType;
import prism.PrismException;

public class StateUpdateFactory
{

	public static StateUpdate getStateUpdate(ModelType model, CachedModelGenerator modelGenerator, StateValueContainer container, int bound, boolean min,
			double epsilon) throws PrismException
	{
		switch (model) {
		case MDP:
			return new StateUpdateMDP(modelGenerator, container, bound, min, epsilon);
			/* TODO: PORT
			case SMG:
				return new StateUpdateSMG(me, container, bound, min, epsilon);
			*/
		}
		throw new PrismException("Unknown model, got " + model);
	}

}
