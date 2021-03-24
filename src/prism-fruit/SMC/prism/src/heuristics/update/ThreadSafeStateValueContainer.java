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

import heuristics.search.StateValue;
import parser.State;

import java.util.concurrent.locks.ReentrantLock;

public final class ThreadSafeStateValueContainer implements StateValueContainer
{
	private final StateValueContainer delegate;
	private final ReentrantLock qLock = new ReentrantLock();
	private final ReentrantLock actionLock = new ReentrantLock();

	public ThreadSafeStateValueContainer()
	{
		this(new DefaultDepthStateValueContainer());
	}

	public ThreadSafeStateValueContainer(StateValueContainer delegate)
	{
		this.delegate = delegate;
	}

	@Override public boolean setQValue(State state, StateValue value, int depth)
	{
		qLock.lock();
		try {
			return delegate.setQValue(state, value, depth);
		} finally {
			qLock.unlock();
		}
	}

	@Override public StateValue getQValue(State state, int depth)
	{
		qLock.lock();
		try {
			return delegate.getQValue(state, depth);
		} finally {
			qLock.unlock();
		}
	}

	@Override public void setCurrentAction(State state, int action)
	{
		actionLock.lock();
		try {
			delegate.setCurrentAction(state, action);
		} finally {
			actionLock.unlock();
		}
	}

	@Override public int getCurrentAction(State state)
	{
		actionLock.lock();
		try {
			return delegate.getCurrentAction(state);
		} finally {
			actionLock.unlock();
		}
	}

	@Override public int getSize()
	{
		return delegate.getSize();
	}

	@Override public int getNumAllValues()
	{
		qLock.lock();
		try {
			return delegate.getNumAllValues();
		} finally {
			qLock.unlock();
		}
	}
}
