#ifndef RBT_HPP
#define RBT_HPP
const bool RBT_RED = true;
const bool RBT_BLACK = false;
template <typename Key, typename T>
class RBNode
{
public:
    Key key;
    T *value;
    size_t n; //sum
    RBNode *left, *right;
    bool color;
    //
    RBNode(Key _key, T *_value, size_t _n, bool _color)
    {
        key = _key;
        value = _value;
        n = _n;
        left = right = nullptr;
        color = _color;
    }
};

template <typename Key, typename T>
class RBT
{
public:
    RBT() { root = nullptr; }
    ~RBT() { destory_(root); }
    int size() { return size_(root); }
    T *get(Key key) { return get_(root, key); }
    void put(Key key, T *value)
    {
        root = put_(root, key, value);
        root->color = RBT_BLACK;
    }
    Key min() { return min_(root)->key; }
    Key max() { return max_(root)->key; }
    void deleteMin()
    {
        if (!isRed_(root->left) && !isRed_(root->right))
            root->color = RBT_RED;
        root = deleteMin_(root);
        if (root)
            root->color = RBT_BLACK;
    }
    void deleteMax_()
    {
        if(!isRed_(root->left) && !isRed_(root->right))
            root->color = RBT_RED;
        root = deleteMax_(root);
        if(root)
            root->color = RBT_BLACK;
    }
    void deleteNode(Key key)
    { 
        if(!isRed_(root->left) && !isRed_(root->right))
            root->color = RBT_RED;
        root = deleteNode_(root, key);
        if(root)
            root->color = RBT_BLACK;
    }

    std::queue<Key> keys() { return keys(min(), max()); }
    std::queue<Key> keys(Key lo, Key hi);

    void print() //for test
    {
        std::cout << "node: ";
        print_(root);
        std::cout << std::endl;
    }

private:
    typedef RBNode<Key, T> Node;
    Node *root;

    void destory_(Node *node)
    {
        if (node)
        {
            destory_(node->left);
            destory_(node->right);
            delete node;
            node = nullptr;
        }
    }

    int size_(Node *node)
    {
        if (node == nullptr)
            return 0;
        else
            return node->n;
    }

    T *get_(Node *node, Key key)
    {
        if (node == nullptr)
            return nullptr;
        if (key < node->key)
            return get_(node->left, key);
        else if (key > node->key)
            return get_(node->right, key);
        else
            return node->value;
    }

    bool isRed_(Node *node)
    {
        if (node == nullptr)
            return false;
        return node->color == RBT_RED;
    }
    Node *rotateLeft_(Node *node)
    {
        Node *t = node->right;
        node->right = t->left;
        t->left = node;
        t->n = node->n;
        node->n = size_(node->left) + size_(node->right) + 1;

        t->color = node->color;
        node->color = RBT_RED;
        return t;
    }
    Node *rotateRight_(Node *node)
    {
        Node *t = node->left;
        node->left = t->right;
        t->right = node;
        t->n = node->n;
        node->n = size_(node->left) + size_(node->right) + 1;

        t->color = node->color;
        node->color = RBT_RED;
        return t;
    }

    void flipColors_(Node *node)
    {
        node->color = RBT_RED;
        node->left->color = RBT_BLACK;
        node->right->color = RBT_BLACK;
    }

    Node *put_(Node *node, Key key, T *value)
    {
        if (node == nullptr)
            return new Node(key, value, 1, RBT_RED);
        if (key < node->key)
            node->left = put_(node->left, key, value);
        else if (key > node->key)
            node->right = put_(node->right, key, value);
        else
            node->value = value;

        if (isRed_(node->right) && !isRed_(node->left))
            node = rotateLeft_(node);
        if (isRed_(node->left) && isRed_(node->left->left))
            node = rotateRight_(node);
        if (isRed_(node->left) && isRed_(node->right))
            flipColors_(node);

        node->n = size_(node->left) + size_(node->right) + 1;
        return node;
    }

    Node *min_(Node *node)
    {
        if (node->left == nullptr)
            return node;
        return min_(node->left);
    }

    Node *max_(Node *node)
    {
        if (node->right == nullptr)
            return node;
        return max_(node->right);
    }

    Node *moveRedLeft_(Node *node)
    {
        flipColors_(node);
        if (isRed_(node->right->left))
        {
            node->right = rotateRight_(node->right);
            node = rotateLeft_(node);
        }
        return node;
    }
    
    Node *moveRedRight_(Node * node)
    {
        flipColors_(node);
        if(!isRed_(node->left->left))
            node = rotateRight_(node);
        return node;
    }

    Node *balance_(Node *node)
    {
        if (isRed_(node->right))
            node = rotateLeft_(node);

        if (isRed_(node->right) && !isRed_(node->left))
            node = rotateLeft_(node);
        if (isRed_(node->left) && isRed_(node->left->left))
            node = rotateRight_(node);
        if (isRed_(node->left) && isRed_(node->right))
            flipColors_(node);

        node->n = size_(node->left) + size_(node->right) + 1;
        return node;
    }

    Node *deleteMin_(Node *node)
    {
        if (node->left == nullptr)
        {
            delete node;
            return nullptr;
        }

        if (!isRed_(node->left) && !isRed_(node->left->left))
            node = moveRedLeft_(node);
        node->left = deleteMin_(node->left);
        return balance_(node);
    }

    Node *deleteMax_(Node *node)
    {
        if(isRed_(node->left))
            node = rotateRight_(node);
        if(node->right == nullptr){
            delete node;
            return nullptr;
        }
        if(!isRed_(node->right) && !isRed_(node->right->left))
            node = moveRedRight_(node);
        node->right = deleteMax_(node->right);
        return balance_(node);
    }

    Node *deleteNode_(Node* node, Key key)
    {
        if(key < node->key)
        {
            if(!isRed_(node->left) && !isRed_(node->left->left))
                node = moveRedLeft_(node);
            node->left = deleteNode_(node->left, key);
        }
        else
        {
            if(isRed_(node->left))
                node = rotateRight_(node);
            if(key == node->key && node->right == nullptr)
            {
                delete node;
                return nullptr;
            }
            if(!isRed_(node->right) && !isRed_(node->right->left))
                node = moveRedRight_(node);
            if(key == node->key)
            {
                node->value = get_(node->right, min_(node->right)->key);
                node->key = min_(node->right)->key;
                node->right = deleteMin_(node->right);
            }
            else
                node->right = deleteNode_(node->right, key);
        }
        return balance_(node);
    }

    void keys_(Node *node, std::queue<Key> &q, Key lo, Key hi)
    {
        if (node == nullptr)
            return;
        if (lo < node->key)
            keys_(node->left, q, lo, hi);
        if (lo <= node->key && hi >= node->key)
            q.push(node->key);
        if (hi > node->key)
            keys_(node->right, q, lo, hi);
    }

    void print_(Node *node) //for test
    {
        if (node)
        {
            print_(node->left);
            std::cout << node->key << " ";
            print_(node->right);
        }
    }
};

template <typename Key, typename T>
inline std::queue<Key> RBT<Key, T>::keys(Key lo, Key hi)
{
    if (hi < lo)
    {
        Key t = lo;
        lo = hi;
        hi = t;
    }

    std::queue<Key> q;
    keys_(root, q, lo, hi);
    return q;
}

#endif // RBT_HPP