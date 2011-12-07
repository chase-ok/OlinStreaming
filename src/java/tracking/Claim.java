
package tracking;

import java.util.*;

public class Claim {
    
    public static class Stake implements Comparable<Stake> {
        
        public Object claimer;
        public double score;
        public int rank;
        
        public Stake(Object claimer, double score, int rank) {
            this.claimer = claimer;
            this.score = score;
            this.rank = rank;
        }
        
        public int compareTo(Stake other) {
            if (rank == other.rank) return Double.compare(score, other.score);
            return other.rank - rank;
        }
        
    }
    
    public Object object;
    public ArrayList<Stake> stakes;
    
    public Claim(Object[] object) {
        this.object = object[0];
        stakes = new ArrayList<Stake>();
    }
    
    public boolean hasStakes() {
        return stakes.size() != 0;
    }
    
    public void stake(Object[] claimer, double score, int rank) {
        stakes.add(new Stake(claimer[0], score, rank));
    }
    
    public void revoke(Object[] claimer) {
        for (int i = 0; i < stakes.size(); i++) {
            if (stakes.get(i).claimer == claimer[0]) {
                stakes.remove(i);
                return;
            }
        }
    }
    
    public Object getWinner() {
        Collections.sort(stakes);
        return stakes.get(0).claimer;
    }
    
}
