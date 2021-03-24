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
import parser.State;
import prism.ModelGenerator;
import prism.PrismException;

public final class ValueInit
{
	private final StateUpdate stateUpdate;

	public ValueInit(StateUpdate su)
	{
		this.stateUpdate = su;
	}

	//Add self loop detection

	public double getLowerBoundValue(ModelGenerator pmg, State s, int depth) throws PrismException
	{
		assert depth >= 0;
		checkDepth(depth);

		if (stateUpdate.isTarget(s))
			return 1;
		if (depth >= stateUpdate.getBound())
			return 0;

		return 0;
	}

	public double getLowerBoundValue(ModelGenerator pmg, State s) throws PrismException
	{
		// TODO This can be omitted / simplify in the MP case.
		if (stateUpdate.isTarget(s))
			return 1;
		if (stateUpdate.isZeroState(s))
			return 0;
		return 0;
	}

	public double getUpperBoundValue(ModelGenerator pmg, State s) throws PrismException
	{
		if (stateUpdate.isTarget(s))
			return 1;
		if (stateUpdate.isZeroState(s))
			return 0;
		return 1;
	}

	public double getUpperBoundValue(ModelGenerator pmg, State s, int depth) throws PrismException
	{
		assert depth >= 0;
		checkDepth(depth);

		if (stateUpdate.isTarget(s))
			return 1;
		if (depth >= stateUpdate.getBound())
			return 0;
		if (stateUpdate.isSelfLoop(s))
			return 0;

		return 1;
	}

	private void checkDepth(int depth) throws PrismException
	{
		if (depth > stateUpdate.getBound() + 1) {
			throw new PrismException("Depth cannot be bigger than bound");
		}
	}

}
