#include "tree-lib.hpp"
#include "read_characters.h"
#include<queue>
#include<algorithm>
#include<string>
#include<cstdint>
#include<vector>
#include<fstream>
#include<list>
#include<array>
#include<unordered_set>
#include<climits>

// Assume in the input tree file there is only contain one newick string.
std::string get_newick(std::string input_tree) {
  std::ifstream file(input_tree);
  
  if (file.fail()) {
    std::cout << "Failed to open the input tree file" << std::endl;
    exit(EXIT_FAILURE);
  }
  std::string line;

  std::getline(file, line);
   
  file.close();
  return line;
}

std::unordered_set<int> intersect_all_sets(const std::vector<std::unordered_set<int>>& sets) {
    
    if (sets.empty()) return {};
    if (sets.size() == 1) return sets[0];
    
    int min_size = 10000000;
    size_t smallest_idx = 0;
    for (size_t i = 1; i < sets.size(); i++) {
        if ((sets[i].size() < min_size) && (sets[i].size() > 0)) {
            smallest_idx = i;
	    min_size = sets[i].size();
        }
    }
    
    std::unordered_set<int> result = sets[smallest_idx];
 	
    for (const int& val : sets[smallest_idx]) {
        bool in_all_sets = true;
        
        for (size_t i = 0; i < sets.size(); i++) {
            if (i == smallest_idx) continue; 
            
            if (sets[i].find(val) == sets[i].end()) {
                in_all_sets = false;
                break;
            }
        }
        
        if (in_all_sets) {
            result.insert(val);
        }
    }
    
    return result;
}

std::string contract_edges(Node *r, std::unordered_set<int> edges){
  
  for (auto j = Traverse::PostOrder(r); j != j.end(); j++) {
    if ((*j)->ID == r->ID){
    	continue;
    }
    else{
	//if node is in the set then contract
	if(edges.find((*j)->ID) != edges.end() && !(*j)->is_leaf()){
		(*j)->contract();
	}
    }
  }	

  return r->newick(false);
}
std::unordered_set<int> one_step(int i, uint8_t** C, unsigned int k, boost::unordered_map<std::string, unsigned int> &label2index, Node* r, size_t max_id) {
	
  std::unordered_set<int> res;
  int *below = new int[max_id];
  int *above = new int[max_id];
  int *character = new int[max_id];
  
  // compute below for every vertex
  for (auto j = Traverse::PostOrder(r); j != j.end(); j++) {
          
    if ((*j)->is_leaf()) {
      std::string leaf_label = (*j)->label;
      //std::cout << "leaf_label: " << leaf_label << std::endl;
      below[(*j)->ID] = C[i][label2index[leaf_label]];
      //std::cout << "below[" << (*j)->ID <<"]" << ": " << below[(*j)->ID] << std::endl;
    } else {
      std::list<Node*> children = (*j)->get_children();
      int missing_num = 0;
      int one_num = 0;
      for (auto child : children) {
	if (below[child->ID] == 2) {
	  missing_num++;
	} else if (below[child->ID] == 1) {
	  one_num++;
	}
      }
      if (missing_num == children.size()) {
	below[(*j)->ID] = 2;
      } else if (one_num >= 1) {
	below[(*j)->ID] = 1;
      } else {
	below[(*j)->ID] = 0;
      }
    }
  }
  //std::cout << "Computed below" << std::endl;
  // compute above for every vertex
  above[r->ID] = -1;
  std::list<Node*> root_children = r->get_children();
  auto it = root_children.begin();
  
  Node* root_left = root_children.front();
  //std::cout << "computed root left" << std::endl;
  std::advance(it, 1);
  Node* root_right = *it;
  above[root_left->ID] = below[root_right->ID];
  //std::cout << "computed root right" << std::endl;
  above[root_right->ID] = below[root_left->ID];
  //std::cout << "computed root left right above" << std::endl;

  for (auto j = Traverse::PreOrder(r); j != j.end(); j++) {

    if ((*j)->ID == (r)->ID || (*j)->ID == root_left->ID || (*j)->ID == root_left->ID) {
      continue;
    }
    
    Node* parent = (*j)->get_parent();
    std::list<Node*> siblings = parent->get_children();
    
    int tot_size_to_check = siblings.size();
    int miss_num = 0;
    int one_num = 0;
    for (Node* sibling : siblings) {

      if (sibling->ID != (*j)->ID) {
	
	if (below[sibling->ID] == 2) {
	  miss_num++;
	} else if (below[sibling->ID] == 1) {
	  one_num++;
	}
      }
    }
    if (above[parent->ID] == 2) {
      miss_num++;
    } else if (above[parent->ID] == 1) {
      one_num++;
    }
    //std::cout << "Compute above:" << (*j)->ID << std::endl;
    if (miss_num == tot_size_to_check) {
      above[(*j)->ID] = 2;
    } else if (one_num >= 1) {
      above[(*j)->ID] = 1;
    } else {
      above[(*j)->ID] = 0;
    }
    //std::cout << "above[" << (*j)->ID << "]: " << above[(*j)->ID] << std::endl;
  }

  //std::cout << "Computed all above" << std::endl;
  // computing character for every vertex
  for (auto j = Traverse::PostOrder(r); j != j.end(); j++) {
    if ((*j)->is_leaf()) {
      character[(*j)->ID] = C[i][label2index[(*j)->label]];
      //std::cout << "computed leaf character[" << (*j)->ID << "]: " << character[(*j)->ID] << std::endl;
    } else {
      int one_num = 0;
      int miss_num = 0;
      std::list<Node*> children = (*j)->get_children();
      //int size_to_check = children.size();

      for (auto child : children) {
	if (below[child->ID] == 2) {
	  miss_num++;
	} else if (below[child->ID] == 1) {
	  one_num++;
	}
      }
      
      if (above[(*j)->ID] == 2) {
	miss_num++;
      } else if (above[(*j)->ID] == 1) {
	one_num++;
      }
      
      if (miss_num >= 2) {
	character[(*j)->ID] = 2;
      } else if (one_num >= 2) {
	character[(*j)->ID] = 1;
      } else {
	character[(*j)->ID] = 0;
      }
      // std::cout << "charactera" << (*j)->ID << "]: " << character[(*j)->ID] << std::endl;
    }
  }
  //std::cout << "computed all character" << std::endl;
  for (auto j = Traverse::PostOrder(r); j != j.end(); j++) {
    if ((*j)->ID == r->ID){
    	continue;
    }
 
    if ((character[(*j)->ID] == 0 && character[((*j)->get_parent())->ID] == 1) || (character[(*j)->ID] == 1 && character[((*j)->get_parent())->ID] == 0)){
    	continue;
    }
    else
    {	    
	res.insert((*j)->ID);
    }
  }

  delete[] character;
  delete[] above;
  delete[] below;
  return res;
}


std::string compute_edge_set(uint8_t** C, unsigned int k,std::string input_tree, boost::unordered_map<std::string, unsigned int> &label2index) {
  std::string newick_str = get_newick(input_tree);
  Tree *T = new Tree(newick_str);
  
  size_t max_id = (T->get_root())->get_max_id();
  std::vector<std::unordered_set<int>> edge_sets(k);
 
  for (int i = 0; i < k; i++) {
    edge_sets.push_back(one_step(i, C, k, label2index, T->get_root(), max_id));
    std::cout.flush();
  }


  std::unordered_set<int> res = intersect_all_sets(edge_sets);

  
  std::string processed_newick = contract_edges(T->get_root(),res);

  delete T;

  for (int i = 0; i < k; i++) {
    delete[] C[i];
  }
  delete[] C;
  
  return processed_newick;
}


int main(int argc, char* argv[])
{

	std::string tree_file = "";
	std::string char_file = "";

	if(argc != 3)
	{
		std::cerr << "Invalid Parameter Count" << std::endl;
		return 1;
	}
	else
	{
		tree_file += argv[1];
		char_file += argv[2];
	}

	std::string newick_file = tree_file + ".contracted";
	std::ofstream output(newick_file);

	if(!output.is_open())
	{
		std::cerr << "File Opening Error" << std::endl;
		return 1;
	}


  	unsigned int k;

  	boost::unordered_map<std::string, unsigned int> label2index;

	std::vector<std::string> labels;

  	boost::unordered_map<std::string, std::string> record;

  	uint8_t** C = read_characters(char_file, k, label2index, labels, record);
  
	std::string processed_newick = compute_edge_set(C, k, tree_file, label2index);


	output << processed_newick << ";";


	return 0;
}
