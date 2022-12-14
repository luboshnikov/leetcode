#include <iostream>
#include <map>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <chrono>
#include <cmath>
using namespace std;

class Node {
public:
    int val;
    Node* next;
    Node* random;

    Node(int _val) {
        val = _val;
        next = nullptr;
        random = nullptr;
    }
};

struct ListNode {
    int val;
    ListNode *next;
    ListNode() : val(0), next(nullptr) {}
    ListNode(int x) : val(x), next(nullptr) {}
    ListNode(int x, ListNode *next) : val(x), next(next) {}
};


struct TreeNode {
    int val;
    TreeNode *left;
    TreeNode *right;
    TreeNode() : val(0), left(nullptr), right(nullptr) {}
    TreeNode(int x) : val(x), left(nullptr), right(nullptr) {}
    TreeNode(int x, TreeNode *left, TreeNode *right) : val(x), left(left), right(right) {}
};

class MyCalendar {
public:

    bool book(int start, int end) {
        auto next = booked.lower_bound(start);
        if (
                (next != booked.end() && next->first < end) ||
                (next != booked.begin() && prev(next)->second > start)){
            return false;
        }
        booked[start] = end;
        return true;
    }

private:
    map<int, int> booked;
};

class Solution {
public:
    static vector<int> twoSum(vector<int>& nums, int target) {
        vector<int> result;
        unordered_map<int, int> negativ;
        for(int i = 0; i < nums.size(); ++i){
            auto it = negativ.find(nums[i]);
            if(it == negativ.end()) {
                negativ[target-nums[i]] = i;
            } else {
                result.push_back(it->second);
                result.push_back(i);
            }
        }
        return result;
    }

    static int lengthOfLongestSubstring(string s) {
        unordered_map<char, int> a;
        int count = 0;
        int j = 0, i = 0;
        while(i < s.length()){
            a[s[i]]++;
            if(a.size() == i-j+1){
                count = max(count, i-j+1);
            } else {
                while(a.size() < i-j+1){
                    a[s[j]]--;
                    if(a[s[j]] == 0){
                        a.erase(s[j]);
                    }
                    j++;
                }
            }
            i++;
        }
        return count;
    }

    static bool isPalindrome(int x) {
        string s = to_string(x);
        int k = 0;
        while(k < s.size()/2){
            if(s[k] != s[s.size() - 1 - k]){
                return false;
            }
            k++;
        }
        return true;
    }

    static int value(char x){
        if(x=='I'){
            return 1;
        }
        if(x=='V'){
            return 5;
        }
        if(x=='X'){
            return 10;
        }
        if(x=='L'){
            return 50;
        }
        if(x=='C'){
            return 100;
        }
        if(x=='D'){
            return 500;
        }
        if(x=='M'){
            return 1000;
        }

        return -1;
    }
    static int romanToInt(string s) {
        int result = value(s.back());
        for(int i = 0; i < s.size()-1; i++){
            if (value(s[i]) < value(s[i+1])) result -= value(s[i]);
            else result += value(s[i]);
        }
        return result;
    }

    static string longestCommonPrefix(vector<string>& strs) {
        string res = strs[0];
        int j;
        for(int i = 1; i < strs.size(); i++){
            j = 0;
            while(j < res.size() and j < strs[i].size()){
                if(res[j] != strs[i][j]){
                    break;
                }
                j++;
            }
            res = res.substr(0, j);
        }
        return res;
    }

    static int numberOfArithmeticSlices(const vector<int>& nums) {
        int count = 0, j = 0;
        for(int i = 2; i < nums.size(); i++){
            if(nums[i] - nums[i-1] == nums[i-1] - nums[i-2])
                count += (i-j+1)-2;
            else
                j = i-1;
        }
        return count;
    }

    static double champagneTower(int poured, int query_row, int query_glass) {
        vector<double> res(101, 0);
        res[0] = poured;
        for(int i = 1; i <= query_row; ++i){
            for(int j = i; j >= 0; --j){
                res[j+1] += res[j] = max(0.0, (res[j] - 1)/2.0);
            }
        }
        return min(1.0, res[query_glass]);
    }


    static ListNode* mergeTwoLists(ListNode* list1, ListNode* list2) {
        ListNode res;
        ListNode* n = &res;
        while(list1 && list2){
            if (list1->val < list2->val){
                n->next = list1;
                list1 = list1->next;
                n = n->next;
            } else {
                n->next = list2;
                list2 = list2->next;
                n = n->next;
            }
        }
        if(list2) {
            n->next = list2;
        } else if (list1){
            n->next = list1;
        }
        return res.next;
    }

    static bool hasCycle(ListNode *head) {
        if(!head)
            return false;
        ListNode* first = head, *second = head;
        while(first->next && first->next->next) {
            second = second->next;
            first = first->next->next;
            if(first == second){
                return true;
            }
        }
        return false;

    }

    static ListNode* deleteDuplicates(ListNode* head) {
        if(!head){
            return head;
        }
        ListNode res;
        ListNode* n = &res;
        while(head && head->next){
            if(head->val != head->next->val){
                n->next = head;
                n = n->next;
                head = head->next;
            }
            else {
                while(head->next && head->val == head->next->val) {
                    head = head->next;
                }
                head = head->next;
            }
        }
        n->next = head;
        return res.next;
    }

    static int removeDuplicates(vector<int>& nums) {
        int count = 1;
        for(int i = 1; i < nums.size(); i++){
            if(nums[i-1] != nums[i]){
                nums[count++] = nums[i];
            }
        }
        return count;
    }

    static int search(vector<int>& nums, int target) {
        int left = 0, right = nums.size() - 1, med;
        while(left <= right) {
            med = (right+left)/2;
            if(nums[med] < target) {
                left = med+1;
            } else if(nums[med] > target) {
                right = med-1;
            } else {
                return med;
            }
        }
        return -1;
    }

    static bool isBadVersion(int version) {return true;}
    static int firstBadVersion(int n) {
        int first = 0, med;
        while(first < n) {
            med = first + (n-first)/2;
            if (isBadVersion(med)) {
                n = med;
            } else {
                first = med + 1;
            }
        }
        return n;
    }

    static int searchInsert(vector<int>& nums, int target) {
        int left = 0, right = nums.size() - 1, med;
        while(left <= right) {
            med = left + (right-left)/2;
            if(nums[med] < target) {
                left = med+1;
            } else if(nums[med] > target) {
                right = med-1;
            } else {
                return med;
            }
        }
        return left;
    }

    static vector<int> sortedSquares(vector<int>& nums) {
        vector<int> res(nums.size());
        int l = 0, r = nums.size() - 1;
        for(int i = res.size() - 1; i >= 0; i--){
            if(abs(nums[l]) < abs(nums[r]))
                res[i] = nums[r]*nums[r--];
            else
                res[i] = nums[l]*nums[l++];
        }
        return res;
    }

    static string minRemoveToMakeValid(string s) {
        vector<int> locat;
        for (int i = 0; i < s.size(); i++){
            if (s[i] == '('){
                locat.push_back(i);
            } else if(s[i] == ')'){
                if(locat.empty() || s[locat.back()] == ')'){
                    locat.push_back(i);
                } else {
                    locat.pop_back();
                }
            }
        }
        int j = 0;
        for(int i = 0, id = 0; i < s.size(); i++){
            if(id < locat.size() && i == locat[id]) {
                id++;
            } else {
                s[j++] = s[i];
            }
        }
        return s.substr(0, j);
    }

    static int threeSumMulti(const vector<int>& arr, int target) {
        vector<size_t> count(101, 0);
        for(const auto& i : arr){
            count[i]++;
        }
        size_t res = 0;
        for(int k = min(100, target); k >= target/3; k--){
            int med = target - k;
            for(int j = min(med, k); j >= med/2; j--){
                int i = med - j;
                if(i < 0 || i > j) continue;
                size_t n = count[i] * count[j] * count[k];
                if(n == 0) continue;
                if (i == k) n = count[i] * (count[i] - 1) * (count[i] - 2) / 6;
                else if (i == j) n = count[i] * (count[i] - 1) / 2 * count[k];
                else if (j == k) n = count[i] * count[j] * (count[j]-1) / 2;
                res += n;
            }
        }
        return (int)(res % 1000000007);
    }

    static int kthSmallest(vector<vector<int>>& matrix, int k) {
        int n = matrix.size(), i, min;
        vector<int> vector1(n, 0);
        while(k != 0){
            min = 1000000001;
            for(int j = 0; j < n; j++){
                if(vector1[j] < n and matrix[j][vector1[j]] < min){
                    i = j;
                    min = matrix[i][vector1[i]];
                }
            }
            k--;
            vector1[i]++;
        }
        return min;
    }

    static int combinationSum4(vector<int>& nums, int target) {
        vector<uint> res(target+1, 0);
        vector<int> sort_nums;
        for(auto i : nums){
            if(i <= target){
                res[i] = 1;
            }
        }
        for(int i = 1; i < target+1; i++){
            if(res[i] != 0){
                sort_nums.push_back(i);
            }
            for(auto j : sort_nums){
                res[i] += res[i-j];
            }
        }
        return res[target];
    }

    static int mirrorReflection(int p, int q) {
        int k = 1;
        uint count = 1;
        int loc = q;
        while(loc != 0 && loc != p){
            loc += q*k;
            count++;
            if(loc > p){
                k = -1;
                loc = p - loc;
            } else if(loc < 0){
                loc *= k;
                k = 1;
            }
        }
        if(count%2 == 0){
            return 2;
        }
        return loc/p;
    }

    static int poorPigs(int buckets, int minutesToDie, int minutesToTest) {
        int process = minutesToTest/minutesToDie;
        if(process == 0){
            return 0;
        }
        process++;
        int pig = 1;
        int res = process;
        while(res < buckets){
            res = pow(process, pig);
            pig++;
        }
        return pig;
    }

    static int countVowelPermutation(int n) {
        unsigned int a = 1, a_old, e = 1, e_old, i = 1, i_old, o = 1, o_old, u = 1, u_old, M = 1e9 + 7;

        for(int k = 1; k < n; k++){
            a_old = a, e_old = e, i_old = i, o_old = o, u_old = u;
            a = e_old;
            e = (a_old + i_old)%M;
            i = (a_old+e_old+o_old+u_old)%M;
            o = (i_old+u_old)%M;
            u = a_old;

        }
        return (int)((a+e+i+o+u)%M);
    }

    static int lengthOfLIS(const vector<int>& nums) {
        map<int, unsigned int> g = {{nums[0], 1}};
        int res = 1;
        for(int i = 1; i < nums.size(); i++){
            auto it = g.lower_bound(nums[i]);
            if(it == g.begin()){
                if(it->first != nums[i]){
                    g[nums[i]] = 1;
                }
            } else {
                int max = 0;
                for(auto j = g.begin(); j != it; j++){
                    if(j->second >= max){
                        max = j->second;
                    }
                }
                g[nums[i]] = ++max;
                if(max > res){
                    res = max;
                }
            }
        }

        return res;
    }

    static int numFactoredBinaryTrees(vector<int>& arr) {
        vector<unsigned long> countTrees(arr.size(), 1);
        sort(arr.begin(), arr.end());

        for(int i = 1; i < arr.size(); i++){
            unsigned long count = 1;
            for(int j = 0; j < i; j++){
                if(arr[i]%arr[j] == 0 && binary_search(arr.begin(), arr.begin()+i, arr[i]/arr[j])){
                    int z = lower_bound(arr.begin(), arr.begin()+i, arr[i]/arr[j]) - arr.begin();
                    count += countTrees[j] * countTrees[z];
                }
                countTrees[i] = count;
            }
        }

        unsigned int res = 0, M = 1e9 + 7;
        for(auto& i : countTrees){
            res = (res + i) % M;
        }
        return (int)res;
    }

    static bool isPalindrome(ListNode* head) {
        ListNode *new_head = head;
        head = head->next;
        new_head->next = nullptr;

        if(head == nullptr){
            return true;
        }

        while(head){
            bool res = false;
            ListNode *comp_old;
            ListNode *comp_new;
            if(head->val == new_head->val){
                comp_old = head->next;
                comp_new = new_head->next;
                while(comp_old and comp_new and comp_old->val == comp_new->val){
                    comp_old = comp_old->next;
                    comp_new = comp_new->next;
                }
                if(comp_old == nullptr and comp_new == nullptr){
                    return true;
                }
            }

            if(head->next and head->next->val == new_head->val) {
                comp_old = head->next->next;
                comp_new = new_head->next;
                while(comp_old and comp_new and comp_old->val == comp_new->val){
                    comp_old = comp_old->next;
                    comp_new = comp_new->next;
                }
                if(comp_old == nullptr and comp_new == nullptr){
                    return true;
                }
            }

            comp_new = head;
            head = head->next;
            comp_new->next = new_head;
            new_head = comp_new;
        }
        return false;
    }

    static bool isPowerOfThree(int n) {
        if(n < 1) {
            return false;
        }
        while(n % 3 == 0){
            n /= 3;
        }
        return n == 1;
    }


    static TreeNode* BST(vector<int>::iterator begin, vector<int>::iterator end){
        if(end - begin == 0){
            return nullptr;
        } else if(end - begin == 1){
            return new TreeNode(*begin);
        }
        auto mid = begin + (end - begin)/2;
        return new TreeNode(*mid, BST(begin, mid), BST(mid+1, end));
    }
    static TreeNode* sortedArrayToBST(vector<int>& nums) {
        return BST(nums.begin(), nums.end());
    }

    static bool reorderedPowerOf2(int n) {
        string s = to_string(n);
        sort(s.begin(), s.end());
        unordered_set<string> powerOfTwo;
        for(int i = 0; i < 32; i++){
            uint k = 1 << i;
            string comp = to_string(k);
            if(comp.size() == s.size()){
                sort(comp.begin(), comp.end());
                powerOfTwo.insert(comp);
            }
        }
        return powerOfTwo.count(s) > 0;
    }

    bool canConstruct(string ransomNote, string magazine) {
        unordered_map<char, uint> chars;
        for(auto& i : magazine){
            chars[i]++;
        }
        for(auto& i : ransomNote){
            if(chars.count(i) == 0 or chars[i]-- < 1){
                return false;
            }
        }
        return true;
    }

    static void mark_land(vector<vector<char>>& grid, size_t x, size_t y){
        if(x < 0 or y < 0 or x >= grid.size() or y >= grid[0].size() or grid[x][y] != '1') return;

        grid[x][y] = '2';

        mark_land(grid, x+1, y);
        mark_land(grid, x, y+1);
        mark_land(grid, x-1, y);
        mark_land(grid, x, y-1);
    }
    static int numIslands(vector<vector<char>>& grid) {
        int count = 0;
        for(size_t i = 0; i < grid.size(); i++){
            for(size_t j = 0; j < grid[0].size(); j++){
                if(grid[i][j] == '1'){
                    count++;
                    mark_land(grid, i, j);
                }
            }
        }
        return count;
    }

    static void rotate(vector<vector<int>>& matrix) {
        int sw;
        size_t sz = matrix.size() - 1;
        for(int i = 0; i < matrix.size()/2; i++){
            for(int j = i; j < sz - i; j++){
                sw = matrix[i][j];
                matrix[i][j] = matrix[sz-j][i];
                matrix[sz-j][i] = matrix[sz-i][sz-j];
                matrix[sz-i][sz-j] = matrix[j][sz-i];
                matrix[j][sz-i] = sw;
            }
        }
    }

    static void checkFlow(vector<vector<int>>& heights, vector<vector<bool>>& land, size_t x, size_t y){
        if(land[x][y]){
            return;
        }
        land[x][y] = true;
        if(x < heights.size()-1 and !land[x+1][y] and heights[x+1][y] >= heights[x][y]){
            checkFlow(heights, land, x+1, y);
        }
        if(y < heights[0].size()-1 and !land[x][y+1] and heights[x][y+1] >= heights[x][y]){
            checkFlow(heights, land, x, y+1);
        }
        if(x > 0 and !land[x-1][y] and heights[x-1][y] >= heights[x][y]){
            checkFlow(heights, land, x-1, y);
        }
        if(y > 0 and !land[x][y-1] and heights[x][y-1] >= heights[x][y]){
            checkFlow(heights, land, x, y-1);
        }
    }

    static vector<vector<int>> pacificAtlantic(vector<vector<int>>& heights) {
        size_t last_row = heights.size(), last_column = heights[0].size();
        vector<vector<bool>> atlantic(last_row, vector<bool>(last_column, false)), pacific(last_row, vector<bool>(last_column, false));
        for(size_t row = 0; row < last_row; row++){
            checkFlow(heights, atlantic, row, last_column-1);
            checkFlow(heights, pacific, row, 0);
        }
        for(size_t column = 0; column < last_column; column++){
            checkFlow(heights, atlantic, last_row - 1, column);
            checkFlow(heights, pacific, 0, column);
        }

        vector<vector<int>> res;
        for(int i = 0; i < last_row; i++){
            for(int j = 0; j < last_column; j++){
                if(atlantic[i][j] and pacific[i][j]){
                    res.push_back({i, j});
                }
            }
        }
        return res;
    }

    static void counting_good_nodes(TreeNode* root, int& count, int max_val){
        if(root){
            if(root->val >= max_val){
                count++;
                max_val = root->val;
            }
            counting_good_nodes(root->left, count, max_val);
            counting_good_nodes(root->right, count, max_val);
        }

    }

    static int goodNodes(TreeNode* root) {
        int result = 0;
        counting_good_nodes(root, result, -10001);
        return result;
    }
};

