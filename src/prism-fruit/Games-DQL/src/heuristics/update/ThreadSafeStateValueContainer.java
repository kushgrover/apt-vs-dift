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

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.locks.ReentrantLock;

import parser.State;
import heuristics.search.StateValue;

public class ThreadSafeStateValueContainer extends StateValueContainer{
	
	protected Map<State,Map<Integer, StateValue>> qValues = new ConcurrentHashMap<State, Map<Integer,StateValue>>();
	protected Map<State,Integer> currentAction = new ConcurrentHashMap<State, Integer>();
	
	private ReentrantLock lock = new ReentrantLock();
	
	
	public int getNumAllValues() {
		lock.lock();
		try {
			return super.getNumAllValues();
		}finally {
			lock.unlock();
		}
	}
	
}
