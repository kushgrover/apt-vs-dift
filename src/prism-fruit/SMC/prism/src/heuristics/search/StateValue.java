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

import prism.PrismUtils;

public final class StateValue
{
	public final double lowerBound;
	public final double upperBound;

	public StateValue(double lowerBound, double upperBound)
	{
		assert lowerBound <= upperBound;
		this.lowerBound = lowerBound;
		this.upperBound = upperBound;
	}

	public String toString()
	{
		return String.format("[%f,%f] (%g)", lowerBound, upperBound, upperBound - lowerBound);
	}

	@Override
	public boolean equals(Object o)
	{
		if (!(o instanceof StateValue)) {
			return false;
		}
		StateValue sv = (StateValue) o;
		return lowerBound == sv.lowerBound && upperBound == sv.upperBound;
	}

	public boolean close(StateValue other, double epsilon)
	{
		return PrismUtils.doublesAreClose(other.lowerBound, lowerBound, epsilon, true)
				&& PrismUtils.doublesAreClose(other.upperBound, upperBound, epsilon, true);
	}
}
