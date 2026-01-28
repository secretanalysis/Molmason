//! MolMason engine_host crate

pub fn hello() -> &'static str {
    "MolMason engine_host"
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_hello() {
        assert!(hello().contains("MolMason"));
    }
}
