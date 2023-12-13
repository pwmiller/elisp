(require 'eieio)

(defclass
 Random ()
 ((w :initform 32 :type number :documentation "Word size")
  (n :initform 624 :type number :documentation "Degree of recurrence")
  (m
   :initform 397
   :type number
   :documentation "Middle word, an offset used in the recurrence relation")
  (r
   :initform 31
   :type number
   :documentation "Separation point of one word, or the number of bits of the lower bitmask")
  (a
   :initform #x9908B0DF
   :type number
   :documentation "Coefficients of the rational normal form twist matrix")
  (u
   :initform 11
   :type number
   :documentation "Additional Mersenne Twister tempering bit shift/mask")
  (d :initform #xFFFFFFFF :type number :documentation "Tempering bitmask d")
  (s :initform 7 :type number :documentation "TGFSR(R) tempering bit shift")
  (b
   :initform #x9D2C5680
   :type number
   :documentation "TGFSR(R) tempering bitmask")
  (t :initform 15 :type number :documentation "Tempering bit shift")
  (c :initform #xEFC60000 :type number :documentation "Tempering bitmask")
  (l :initform 18 :type number :documentation "Additional tempering bit shift")
  (f
   :initform 1812433253
   :type number
   :documentation "Initialization multiplier")
  (MT
   :initform (make-vector 624 0)
   :type vector
   :documentation "The array for the state vector")
  (index
   :initform 625
   :type number
   :documentation "Current index in the state vector")
  (lower_mask :initform #x7FFFFFFF :type number :documentation "Lower bitmask")
  (upper_mask :initform #x80000000 :type number :documentation "Upper bitmask")
  (default_seed :initform 5489 :type number :documentation "Default seed")
  (c_seed
   :initarg
   :c_seed
   :type number
   :documentation "Current seed for he generator"))
 "A class representing a Mersenne Twister pseudorandom number generator.")

(cl-defmethod initialize-instance :after
  ((obj Random) &rest args)
  "Initialize the Random object."
  ;; Set constants like w, n, m, r, a, u, d, s, b, t, c, l, f
  (oset obj w 32)
  (oset obj n 624)
  (oset obj m 397)
  (oset obj r 31)
  (oset obj a #x9908B0DF)
  (oset obj u 11)
  (oset obj d #xFFFFFFFF)
  (oset obj s 7)
  (oset obj b #x9D2C5680)
  (oset obj t 15)
  (oset obj c #xEFC60000)
  (oset obj l 18)
  (oset obj f 1812433253)

  ;; Initialize the state array MT
  (oset obj MT (make-vector (oref obj n) 0))
  (oset obj index (1+ (oref obj n)))

  ;; Initialize masks
  (oset obj lower_mask #x7FFFFFFF)
  (oset obj upper_mask #x80000000)

  ;; Initialize the seed

  (when (not (alist-get :c_seed args))
    (warn "Mersenne Twister seed not provided, using default value.")
    (oset obj c_seed (oref obj default_seed)))
  (seed obj (oref obj c_seed))

  (condition-case nil
      (cl-call-next-method)
    (error nil)))

(cl-defmethod seed ((obj Random) num)
  "Initialize the generator from a seed."
  ;; Set the first element of the MT array to the seed
  (aset (oref obj MT) 0 num)
  (oset obj index (oref obj n))
  ;; Initialize the rest of the elements of the MT array
  (cl-loop
   for i from 1 below (oref obj n) do
   (let* ((prev (aref (oref obj MT) (1- i)))
          (temp
           (+ (* (oref obj f) (logxor prev (lsh prev (- (oref obj w) 2)))) i)))
     (aset (oref obj MT) i (logand temp #xFFFFFFFF)))))


(cl-defmethod twist ((obj Random))
  "Generate the next n values from the series x_i."
  (cl-loop
   for i from 0 below (oref obj n) do
   (let* ((x
           (+ (logand (aref (oref obj MT) i) (oref obj upper_mask))
              (logand
               (aref (oref obj MT) (mod (1+ i) (oref obj n)))
               (oref obj lower_mask))))
          (xA (lsh x -1)))
     ;; Apply the twist transformation
     (when (/= (mod x 2) 0)
       (setq xA (logxor xA (oref obj a))))
     (aset
      (oref obj MT) i
      (logxor (aref (oref obj MT) (mod (+ i (oref obj m)) (oref obj n))) xA)))
   (oset obj index 0)))

(cl-defmethod extract-number ((obj Random))
  "Extract a tempered value based on MT[index], calling twist() every n numbers."
  ;; Call twist if necessary
  (when (>= (oref obj index) (oref obj n))
    (twist obj))

  ;; Extract the value
  (let ((y (aref (oref obj MT) (oref obj index))))
    ;; Apply tempering
    (setq y (logxor y (lsh y (- (oref obj u)))))
    (setq y (logxor y (logand (lsh y (oref obj s)) (oref obj b))))
    (setq y (logxor y (logand (lsh y (oref obj t)) (oref obj c))))
    (setq y (logxor y (lsh y (- (oref obj l)))))

    ;; Increment index and return tempered value
    (oset obj index (1+ (oref obj index)))
    (logand y #xFFFFFFFF)))

(cl-defmethod mt-random ((obj Random))
  "Return a random number in the range [0, 1)."
  (/ (extract-number obj) (expt 2.0 32)))

(cl-defmethod randint ((obj Random) a b)
  "Return a random integer in the range [a, b)."
  (let ((range (- b a))
        (random-value (mt-random obj)))
    (+ a (floor (* random-value range)))))

(cl-defmethod shuffle ((obj Random) list)
  "Shuffle the given list."
  (let ((new-list (copy-sequence list))
        (i (1- (length list))))
    (while (> i 0)
      (let* ((j (randint obj 0 (1+ i)))
             (temp (nth i new-list)))
        (setf (nth i new-list) (nth j new-list))
        (setf (nth j new-list) temp))
      (setq i (1- i)))
    new-list))
