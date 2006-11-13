#ifndef __TreeUtilities_h__
#define __TreeUtilities_h__

#include <stack>

template<class T, class S>
void findAndErase( T& container, S& item )
{
	T new_container;
	for( typename T::iterator t_iter = container.begin(); t_iter != container.end(); t_iter++ )
		if( *t_iter != item )
			new_container.push_back( *t_iter );
	container = new_container;
};

/**
 * Depth first search to check whether a subtree contains a given node
 */
template<class Tree>
bool containsNode( Tree& t, node_id_t subtree_nodeI, node_id_t query_nodeI )
{
	std::stack< node_id_t > node_stack;
	node_stack.push( subtree_nodeI );
	while( node_stack.size() > 0 )
	{
		node_id_t cur_node = node_stack.top();
		node_stack.pop();
		if( cur_node == query_nodeI )
			return true;
		if( t[cur_node].children.size() > 0 )
		{
			for( size_t childI = 0; childI < t[cur_node].children.size(); childI++ )
				node_stack.push( t[cur_node].children[childI] );
		}
	}
	return false;
}


/** place a root on the branch with endpoints root_left and root_right
 */
template<class Tree>
void rerootTree( Tree& t, node_id_t new_root )
{
	// new root must be an internal node
	if( t[new_root].children.size() == 0 )
		throw "Can't root on a leaf node";
	if( new_root == t.root )
		return;	// idiot caller didn't realize it's already rooted here

	// change the old root node to an internal node
	uint childI = 0;
	for( ; childI < t[t.root].children.size(); childI++ ){
		if( containsNode( t, t[t.root].children[childI], new_root ) )
		{
			t[t.root].parents.push_back( t[t.root].children[childI] );
			findAndErase( t[t.root].children, t[t.root].children[childI] );
			break;
		}
	}
	// shake the tree out on the new root node
	t.root = new_root;
	t[t.root].children.insert( t[t.root].children.end(), t[t.root].parents.begin(), t[t.root].parents.end() );
	t[t.root].parents.clear();

	std::stack<node_id_t> node_stack;
	node_stack.push(t.root);
	while( node_stack.size() > 0 )
	{
		// delete the current node from all of its child nodes lists 
		// and insert it as a parent
		// make all other nodes reference by the child grandchildren
		// recurse on each child
		node_id_t cur_node = node_stack.top();
		node_stack.pop();
		for( uint childI = 0; childI < t[cur_node].children.size(); childI++ )
		{
			findAndErase( t[t[cur_node].children[childI]].children, cur_node );
			findAndErase( t[t[cur_node].children[childI]].parents, cur_node );
			t[t[cur_node].children[childI]].children.insert( t[t[cur_node].children[childI]].children.end(), t[t[cur_node].children[childI]].parents.begin(), t[t[cur_node].children[childI]].parents.end() );
			t[t[cur_node].children[childI]].parents.clear();
			t[t[cur_node].children[childI]].parents.push_back(cur_node);
			node_stack.push(t[cur_node].children[childI]);
		}
	}
}


#endif // __TreeUtilities_h__
