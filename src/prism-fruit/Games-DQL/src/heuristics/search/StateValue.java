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

public class StateValue {
	
	private double lowerBound = -1;
	private double upperBound = -1;
	
	public StateValue(double l, double u) {
		lowerBound = l;
		upperBound = u;
	}
	
	public void setLowerBound(double l) {
		lowerBound = l;
	}
	
	public void setUpperBound(double u) {
		upperBound = u;
	}
	
	public double getLowerBound() {
		return lowerBound;
	}
	
	public double getUpperBound() {
		return upperBound;
	}
	
	public String toString() {
		return "Lower bound " + lowerBound + " Upper bound " + upperBound + " Diff " + (upperBound - lowerBound);
	}
	
	public boolean equals(Object o) {
		if(o == null) throw new NullPointerException();
		if(!(o instanceof StateValue)) return false;
		StateValue sv = (StateValue)o;
		boolean equals = true;
		equals = equals && PrismUtils.doublesAreClose(sv.lowerBound, lowerBound, 10e-8, true);
		equals = equals && PrismUtils.doublesAreClose(sv.upperBound, upperBound, 10e-8, true);
		return equals;
		
		
	}
	
}
